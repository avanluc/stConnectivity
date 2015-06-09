#pragma once

#include "prefixSumAsm.cu"
#include "GlobalSync.cu"

__device__ int GlobalCounter;

extern __shared__ unsigned char SMem[];

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



__device__ __forceinline__ void FrontierReserve_Warp(int* GlobalCounter, int founds, int& n, int &totalWarp, int& globalOffset) {
		n = founds;
		totalWarp = warpExclusiveScan<32>(n);
		int oldCounter;
		if (LaneID() == 0)	// && totale != 0)
			oldCounter = atomicAdd(GlobalCounter, totalWarp);

		globalOffset = __shfl(oldCounter, 0);
}



__device__ __forceinline__ void Write(int* devFrontier, int* GlobalCounter, int* Queue, int founds) {
		
		int n, total, globalOffset;
		FrontierReserve_Warp(GlobalCounter, founds, n, total, globalOffset);

		const int pos = globalOffset + n;
		for (int i = 0; i < founds; i++){
			devFrontier[pos + i] = Queue[i];
			//printf("thread %d writing %d at position %d\n", Tid, Queue[i], pos + i);
		}
}



__device__ __forceinline__ void swapDev(int*& A, int*& B) {
	int* temp = A;	// frontiers swap
	A = B;
	B = temp;
}



__global__ void BFS_BlockKernel (	const int* __restrict__	devNode,
									const int* __restrict__	devEdge,
									int2* __restrict__	devDistance,
									int* __restrict__	devSource,
									const int Nsources) {
	int Queue[REG_QUEUE];
	int level = 0;
	int FrontierSize = Nsources;

	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* SMemF2 = (int*) &SMem[F2_OFFSET];

	if (Tid < FrontierSize)
		SMemF1[Tid] = devSource[Tid]; 

	while (FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT) {

		int founds = 0;
		for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE) {
			const int index = SMemF1[t];
			printf("Visiting node %d\n", index );
			const int start = devNode[index];
			int end = devNode[index + 1];

			for (int k = start; k < end; k++) {
				const int dest = devEdge[k];

				if (devDistance[dest].x == INT_MAX) {
					int2 temp = make_int2(level, 0);
					atomicStore(&devDistance[dest], temp);
					//devDistance[dest].x = level;
					Queue[founds++] = dest;
				}
			}
		}

		int WarpPos, n, total;
		Write(SMemF2, &GlobalCounter, Queue, founds); 	//  Util/GlobalWrite.cu

		swapDev(SMemF1, SMemF2);
		level++;
		__syncthreads();
		FrontierSize = GlobalCounter;
	}
}
