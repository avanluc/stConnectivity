#pragma once

#include <assert.h>
#include "prefixSumAsm.cu"
#include "GlobalSync.cu"
#include "definition.cuh"

__device__ int GlobalCounter = 0;
__device__ int globalMax = 0;

extern __shared__ unsigned char SMem[];


 /*
* Assert for CUDA functions
*/
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}



/*
* Self made atomic function to store a int2 value 
*/
__device__ __forceinline__ void atomicStore(int2* address, int2 val){
    unsigned long long* addr_as_ull = (unsigned long long*)address;
    unsigned long long  old = *addr_as_ull;
    unsigned long long  assumed;
    assumed = old;
    atomicCAS(addr_as_ull, assumed, *(unsigned long long*)&val);
    return;
}



/*__device__ __forceinline__ void FrontierReserve_Warp(int* GlobalCounter, int founds, int& n, int &totalWarp, int& globalOffset) {
		n = founds;
		totalWarp = warpExclusiveScan<32>(n);
		int oldCounter;
		if (LaneID() == 0)
			oldCounter = atomicAdd(GlobalCounter, totalWarp);

		globalOffset = __shfl(oldCounter, 0);
}*/



__device__ __forceinline__ void FrontierReserve_Block(int* Front_size, int founds, int& n, int &totalBlock, int& globalOffset) {
	int* SM = (int*) &SMem[TEMP_POS];
	n = founds;
	const int warpId = WarpID();
	SM[warpId] = warpExclusiveScan<32>(n);

	__syncthreads();
	if (Tid < BLOCK_SIZE / 32) {
		int sum = SM[Tid];
		const int total = warpExclusiveScan<BLOCK_SIZE / 32>(sum);

		if (Tid == 0) {
			SM[32] = total;
			SM[33] = atomicAdd(Front_size, total);
		}
		SM[Tid] = sum;
	}
	__syncthreads();

	n += SM[warpId];
	totalBlock = SM[32];
	globalOffset = SM[33];
}



__device__ __forceinline__ void Write(int* devFrontier, int* Front_size, int* Queue, int founds) {
		
		int n, total, globalOffset;
		FrontierReserve_Block(Front_size, founds, n, total, globalOffset);

		const int pos = globalOffset + n;
		for (int i = 0; i < founds; i++){
			cudaAssert((pos + i) < BLOCK_FRONTIER_LIMIT, (pos + i));
			if((pos + i) < BLOCK_FRONTIER_LIMIT)
				devFrontier[pos + i] = Queue[i];
		}
}



__device__ __forceinline__ void swapDev(int*& A, int*& B) {
	int* temp = A;	// frontiers swap
	A = B;
	B = temp;
}



__global__ void BFS_BlockKernel (	const int* __restrict__	devNode,
									const int* __restrict__	devEdge,
									const int* __restrict__	devSource,
									int2* __restrict__	devDistance,
									bool* __restrict__ Matrix,
									const int Nsources) {
	int Queue[REG_QUEUE];
	int level = 0;
	//int FrontierSize = Nsources;
	int FrontierSize = 1;

	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	//int* SMemF2 = (int*) &SMem[F2_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];

	if (Tid < FrontierSize)
		SMemF1[Tid] = devSource[blockIdx.x]; 
		//SMemF1[Tid] = devSource[Tid]; 

	while (FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT )
	{

		int founds = 0;
		for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE){
			const int index = SMemF1[t];
			const int start = devNode[index];
			const int2 current = devDistance[index];
			int end = devNode[index + 1];	

			for (int k = start; k < end; k++) {
				const int dest = devEdge[k];
				const int2 destination = devDistance[dest];	

				/* add dest to the queue */
				if(ATOMIC)
				{	
					int old = atomicCAS(&devDistance[dest].x, INT_MAX, level);
					if ( old == INT_MAX) {	
						devDistance[dest].x = level;
						devDistance[dest].y = current.y;
						Queue[founds++] = dest;
					}
					/* update adj matrix */
					else if (destination.y != current.y && destination.y < Nsources){	
						Matrix[ (current.y * Nsources) + destination.y ] = true;	
						Matrix[ (destination.y * Nsources) + current.y ] = true;	
					}
				}
				else if(BFS)
				{
					if (destination.x == INT_MAX) {	
						devDistance[dest].x = level;
						devDistance[dest].y = current.y;
						Queue[founds++] = dest;
					}
				}
				else
				{
					if (destination.x == INT_MAX) {	
						devDistance[dest].x = level;
						devDistance[dest].y = current.y;
						Queue[founds++] = dest;
					}
					else if (destination.y != current.y && destination.y < Nsources){	
						Matrix[ (current.y * Nsources) + destination.y ] = true;	
						Matrix[ (destination.y * Nsources) + current.y ] = true;	
					}
				}
			}
		}
		
		//int WarpPos, n, total;
		Write(SMemF1, &F2SizePtr[0], Queue, founds);

		//swapDev(SMemF1, SMemF2);
		level++;

		FrontierSize = F2SizePtr[0];
		//if(Tid == 0)
		//	atomicAdd(&GlobalCounter, FrontierSize);

		__syncthreads();
		F2SizePtr[0] = 0;
		//if(Tid == 0 && FrontierSize > globalMax);
			//atomicCAS(&globalMax, globalMax, FrontierSize);
	}
}



// __global__ void BFS_BlockKernel1 (	const int* __restrict__	devNode,
// 									const int* __restrict__	devEdge,
// 									const int* __restrict__	devSource,
// 									int* __restrict__	devDistance,
// 									int* __restrict__	devColor,
// 									bool* __restrict__ Matrix,
// 									const int Nsources) {
// 	int Queue[REG_QUEUE];
// 	int level = 0;
// 	int FrontierSize = Nsources;

// 	int* SMemF1 = (int*) &SMem[F1_OFFSET];
// 	int* SMemF2 = (int*) &SMem[F2_OFFSET];
// 	int* F2SizePtr = (int*) &SMem[F2Size_POS];

// 	if (Tid < FrontierSize)
// 		SMemF1[Tid] = devSource[Tid]; 

// 	while (FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT ) {

// 		int founds = 0;
// 		for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE) {
// 			const int index = SMemF1[t];
// 			const int start = devNode[index];
// 			const int currCol = devColor[index];
// 			int end = devNode[index + 1];

// 			for (int k = start; k < end; k++) {
// 				const int dest = devEdge[k];
// 				const int destCol = devColor[dest];			

// 				/* add dest to the queue */
// 				//int old = atomicCAS(&devDistance[dest], INT_MAX, level);
// 				//if (old == INT_MAX) {	
// 				if (devDistance[dest] == INT_MAX) {	
// 					devDistance[dest] = level;
// 					devColor[dest] = currCol;
// 					Queue[founds++] = dest;
// 				}
// 				/* update adj matrix */
// 				else if (destCol != currCol && destCol < Nsources){	
// 					Matrix[ currCol*Nsources + destCol ] = true;	
// 					Matrix[ destCol*Nsources + currCol ] = true;	
// 				}
// 			}
// 		}
		
// 		int WarpPos, n, total;
// 		Write(SMemF2, &F2SizePtr[0], Queue, founds);

// 		swapDev(SMemF1, SMemF2);
// 		level++;
// 		//__syncthreads();

// 		FrontierSize = F2SizePtr[0];

// 		__syncthreads();
// 		F2SizePtr[0] = 0;
// 	}
// }