#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "CacheFunc.cu"


__global__ void initBitMask(int* BitMask, const int N, const int MaskSize)
{
	if(GTid < (MaskSize - N))
			markWrite<cub::LOAD_CS, cub::STORE_CS, int>(BitMask, N + GTid);
}
	


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
									int* __restrict__ BitMask, 
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
								//BitMask[dest] = 1;									
								markWrite<cub::LOAD_CS, cub::STORE_CS, int>(BitMask, dest);
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
		}
		FrontierSize = 1; 
		level = 0; 
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
	//int Queue[8];

	//int Limit = N % 8 == 0 ? N/8 : (N/8) + 1;

	while( FrontierSize )
	{
		/*if(GTid * 4 < Limit)
		{*/
			int founds = 0;

			/*int* devArray = BitMask + GTid * 4;
			reinterpret_cast<int4*>(Queue)[0] = reinterpret_cast<int4*>( devArray )[0];
			reinterpret_cast<int4*>(Queue)[4] = __ldg(&reinterpret_cast<int4*>( devArray )[Limit*4]);
			*/
			
		for (int index = GTid; index < N; index += MAX_CONCURR_TH)
		{
			//int index = i < 4 ? (GTid*4)+i : ((GTid + Limit)*4)+i;
			
			if(devDistance[index].x == INT_MAX)
			{
				const int start = devNode[index];
				const int end = devNode[index + 1];
				for (int k = start; k < end; k++)
				{	
					const int dest = devEdge[k];
					const int2 destination = devDistance[dest];

					if(markAccess<cub::LOAD_CS, int>(BitMask, dest))
					{
						devDistance[index].x = destination.x + 1;
						devDistance[index].y = destination.y;
						markWrite<cub::LOAD_CS, cub::STORE_CS, int>(BitMask, dest);
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


/*
* KERNEL FUNCTION FOR STCONNECTIVITY
*/
__global__ void STCONN_BlockKernel1 (const int* __restrict__ devNode,
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