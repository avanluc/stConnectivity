#pragma once

#include <../../cub/cub.cuh>

__device__ unsigned int GSync[MAX_CONCURR_BL(BLOCK_SIZE)];
__device__ bool devNextLevel[4];

__global__ void GReset() {
	int bx = blockIdx.x;
    int tx = threadIdx.x;
    int id = tx + (bx*BLOCK_SIZE);
	if (id < MAX_CONCURR_BL(BLOCK_SIZE))
		GSync[id] = 0;
}

/*
*	Glocal Syncronization function
*/
__device__  __forceinline__ void GlobalSync() {
	volatile unsigned *VolatilePtr = GSync;
	__syncthreads();
	
	if (blockIdx.x == 0) {
		if (threadIdx.x == 0){
			VolatilePtr[blockIdx.x] = 1;
		}
		//__syncthreads();

		if (threadIdx.x < MAX_CONCURR_BL(BLOCK_SIZE))
			while ( cub::ThreadLoad<cub::LOAD_CG>(GSync + threadIdx.x) == 0 );

		__syncthreads();

		if (threadIdx.x < MAX_CONCURR_BL(BLOCK_SIZE)){
			VolatilePtr[threadIdx.x] = 0;
		}
	}
	else {
		if (threadIdx.x == 0) {
			VolatilePtr[blockIdx.x] = 1;
			while (cub::ThreadLoad<cub::LOAD_CG>(GSync + blockIdx.x) == 1);
		}
		__syncthreads();
	}
}