#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"


__device__ int GlobalCounter = 0;
__device__ int exitFlag = 0;
__device__ int connected = 0;
__device__ int GQueue[REG_QUEUE];


/*
* BFS KERNEL
*/
__global__ void BFS_BlockKernel (	const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									const int* __restrict__ devSource,
									int2* __restrict__ devDistance,
									bool* __restrict__ Matrix,
									const int Nsources) 
{
	int Queue[REG_QUEUE];
	int FrontierSize = 1;
	int level = 0;

	#if(!ATOMIC)
		exitFlag = 0;
	#endif

	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];
	int* Block_Exit = (int*) &SMem[EXIT_FLAG];

	for(int j = Bid; j < Nsources; j += gridDim.x)
	{
		Block_Exit[0] = 0;
		if (Tid < FrontierSize)
			SMemF1[Tid] = devSource[j];
	
		while ( FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT )
		{
			if(!ATOMIC && exitFlag)
				break;
			if(Block_Exit[0]){
				//printf("Block %d exit\n", Bid );
				break;
			}
	
			int founds = 0; int counter = 0;
			for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE)
			{
				const int index = SMemF1[t];
				const int start = devNode[index];
				const int2 current = devDistance[index];
				int end = devNode[index + 1];	
	
				for (int k = start; k < end; k++)
				{
					const int dest = devEdge[k];
					const int2 destination = devDistance[dest];	
					
					#if ATOMIC
						if(founds < REG_QUEUE)
						{
							counter++;
							if ( atomicCAS(&devDistance[dest].x, INT_MAX, level) == INT_MAX ) {
								devDistance[dest].y = current.y;
								Queue[founds++] = dest;
							}
						}
					#else
						if (destination.x == INT_MAX){	
							devDistance[dest].x = level;
							devDistance[dest].y = current.y;
							Queue[founds++] = dest;
						}
					#endif
				}
			}
			
			Write(SMemF1, &F2SizePtr[0], Queue, Block_Exit, founds);
	
			level++;
			FrontierSize = F2SizePtr[0];

			#if ATOMIC
				atomicAdd(&GlobalCounter, counter);
				/*if (Tid == 0)
					atomicAdd(&GlobalCounter, FrontierSize);*/
			#endif
	
			__syncthreads();
			F2SizePtr[0] = 0;
		}
		level = 0;
		FrontierSize = 1;
	}
}
