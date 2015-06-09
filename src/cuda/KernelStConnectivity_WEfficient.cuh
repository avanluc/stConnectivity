#pragma once

#include "prefixSumAsm.cu"
#include "GlobalSync.cu"
#include <assert.h>

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

__device__ __forceinline__ void FrontierReserve_Block(int* GlobalCounter, int founds, int& n, int &totalBlock, int& globalOffset) {
	int* SM = (int*) SMem;
	n = founds;
	const int warpId = WarpID();
	SM[warpId] = warpExclusiveScan<32>(n);

	__syncthreads();
	if (Tid < BLOCK_SIZE / 32) {
		int sum = SM[Tid];
		const int total = warpExclusiveScan<BLOCK_SIZE / 32>(sum);

		if (Tid == 0) {
			SM[32] = total;
			SM[33] = atomicAdd(GlobalCounter, total);
		}
		SM[Tid] = sum;
	}
	__syncthreads();

	n += SM[warpId];
	totalBlock = SM[32];
	globalOffset = SM[33];
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
		FrontierReserve_Block(GlobalCounter, founds, n, total, globalOffset);

		const int pos = globalOffset + n;
		for (int i = 0; i < founds; i++){
			//printf("\t\twriting %d at position %d\n", Queue[i], pos+i);
			assert((pos + i) < BLOCK_FRONTIER_LIMIT);
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
			const int start = devNode[index];
			int end = devNode[index + 1];

			printf("Thread %d visiting node %d\n", Tid, index );
			for (int k = start; k < end; k++) {
				const int dest = devEdge[k];
				//printf("   Thread %d vede nodo %d\n", Tid, dest);

				if (devDistance[dest].x == INT_MAX) {
					atomicCAS(&devDistance[dest].x, INT_MAX, level);
					if(devDistance[dest].x == level)
					//devDistance[dest].x = level;
						Queue[founds++] = dest;
						printf("\tThread %d aggiunge il nodo %d\n", Tid, dest );
				}
			}
		}

		int WarpPos, n, total;
		Write(SMemF2, &GlobalCounter, Queue, founds); 	//  Util/GlobalWrite.cu

		swapDev(SMemF1, SMemF2);
		level++;
		__syncthreads();
		FrontierSize = GlobalCounter;
		if (Tid == 0)
			printf("Livello %d: and FrontierSize = %d < %d\n", level, FrontierSize, BLOCK_FRONTIER_LIMIT);
		__syncthreads();
		GlobalCounter = 0;
	}
}
