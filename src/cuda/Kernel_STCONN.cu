#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"


/*
* UPDATE BITMASK WITH NEW VISITED NODES
*/
__device__ __forceinline__ void UpdateBitMask(int* BitMask,	int* Frontier, const int FrontierSize)
{
	int limit = FrontierSize < BLOCK_FRONTIER_LIMIT ? FrontierSize : BLOCK_FRONTIER_LIMIT;
	for (int t = Tid; t < limit; t += BLOCK_SIZE)
		BitMask[Frontier[t]] = 1;
}



/*
* KERNEL FUNCTION FOR STCONNECTIVITY
*/
__global__ void STCONN_BlockKernel (const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									const int* __restrict__ devSource,
									int2* __restrict__ devDistance,
									bool* __restrict__ Matrix,
									int*  __restrict__ BitMask, 
									const int Nsources)
{
	int Queue[REG_QUEUE];
	int FrontierSize = 1;
	int level = 0;
	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];

	for(int j = Bid; j < Nsources; j += gridDim.x)
	{
		if (Tid < FrontierSize)
			SMemF1[Tid] = devSource[j];
	
		while ( FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT )
		{
			int founds = 0; int counter = 0;
			for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE)
			{
				const int index = SMemF1[t];
				const int start = devNode[index];
				const int2 current = devDistance[index];
				const int end = devNode[index + 1];	
	
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
								BitMask[dest] = 1;
							}
							else if (destination.y != current.y && destination.y < Nsources){	
								Matrix[ (current.y     * Nsources) + destination.y ] = true;	
								Matrix[ (destination.y * Nsources) + current.y 	   ] = true;	
							}
						}
					#else
						if (destination.x == INT_MAX) {	
							devDistance[dest].x = level;
							devDistance[dest].y = current.y;
							Queue[founds++] = dest;
						}
						else if (destination.y != current.y && destination.y < Nsources){	
							Matrix[ (current.y     * Nsources) + destination.y ] = true;	
							Matrix[ (destination.y * Nsources) + current.y 	   ] = true;	
						}
					#endif
				}
			}

			Write(SMemF1, &F2SizePtr[0], Queue, founds);
	
			level++;
			FrontierSize = F2SizePtr[0];

			#if ATOMIC
				//atomicAdd(&GlobalCounter, counter);
				GlobalWrite(counter, &GlobalCounter);
			#endif

			__syncthreads();
			F2SizePtr[0] = 0;
			//UpdateBitMask(BitMask, SMemF1, FrontierSize);
		}
		//UpdateBitMask(BitMask, SMemF1, FrontierSize);

		level = 0; FrontierSize = 1;
	}

}



/*
* KERNEL FOR BOTTOM-UP VISIT
*/
__global__ void Bottom_Up_Kernel(	const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									int2* __restrict__ devDistance,
									int*  __restrict__ BitMask,
									const int N)
{
	int FrontierSize = 1;

	while( FrontierSize )
	{
		int founds = 0;
		for (int index = GTid; index < N; index += MAX_CONCURR_TH)
		{
			if(devDistance[index].x == INT_MAX)
			{
				const int start = devNode[index];
				const int end = devNode[index + 1];
				for (int k = start; k < end; k++)
				{	
					const int dest = devEdge[k];
					const int2 destination = devDistance[dest];

					if(BitMask[dest] == 1)
					{
						devDistance[index].x = destination.x + 1;
						devDistance[index].y = destination.y;
						BitMask[index] = 1;
						founds++;
						//BottomUp_FrontSize = 1;
						break;
					}
				}
			}
		}
		
		GlobalWrite( founds, &BottomUp_FrontSize);
		
		GlobalSync();
		FrontierSize = BottomUp_FrontSize;
		
		//atomicAdd(&GlobalCounter, founds);
		
		GlobalSync();
		
		BottomUp_FrontSize = 0;
	}
}
