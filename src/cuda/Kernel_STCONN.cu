#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "CacheFunc.cu"
#include <curand_kernel.h>


__global__ void setup_curand(curandState* state)
{
	if(Tid == 0)
		curand_init(Bid, 0, 0, &state[Bid]);
}



__global__ void CheckVisit(	bool* __restrict__ BitMask,
							const int N)
{
	int founds = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH)
		if(BitMask[index] == 0)
			founds++;
	atomicAdd(&VisitResult, founds);
}



__device__ __forceinline__ int visit(	const int* __restrict__ devNode,
										const int* __restrict__ devEdge,
										bool* __restrict__ BitMask,
										const int index)
{
	const int start = devNode[index];
	const int end = devNode[index + 1];
	for (int k = start; k < end; k++){
		if(BitMask[devEdge[k]] == 1){
			return 1;
		}
	}
	return 0;
}



/*
* KERNEL FUNCTION FOR STCONNECTIVITY
*/
__global__ void TopDown_Kernel (const int* __restrict__ devNode,
								const int* __restrict__ devEdge,
								const int* __restrict__ devSource,
								int2* __restrict__ devDistance,
								bool* __restrict__ Matrix,
								bool* __restrict__ BitMask, 
								const double Treshold,
								const int E,
								const int N,
								curandState* globalState)
{
	int Queue[REG_QUEUE];
	int FrontierSize = 1;
	int level = 0;
	uint src  = 0;
	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];

	if(Bid < 2 && Tid == 0){
		color = 2;
		SMemF1[Tid] = devSource[Bid];
		BitMask[SMemF1[Tid]] = 1;
		devDistance[SMemF1[Tid]] = make_int2(0, Bid);
	}

	while(__int2double_rn(GlobalCounter) / __int2double_rn(E) < Treshold)
	{
		if(Tid == 0 && src == 0)
		{
			do{
				src = ((uint)curand(&globalState[Bid])) % N;
			}while(BitMask[src]);
			devDistance[src] = make_int2(0, atomicAdd(&color, 1));
			BitMask[src] = 1;
			SMemF1[Tid] = src;
		}
			
		while ( FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT )
		{
			int founds = 0; int counter = 0;
			for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE)
			{
				const int index = SMemF1[t];
				const int start = devNode[index];
				const int   end = devNode[index + 1];	
				const int2 current = devDistance[index];
	
				for (int k = start; k < end; k++)
				{
					const int dest = devEdge[k];
					const int2 destination = devDistance[dest];	
						
					if(founds < REG_QUEUE)
					{
						counter++;
						if ( atomicCAS(&devDistance[dest].x, INT_MAX, level) == INT_MAX ) {
							devDistance[dest].y = current.y;
							Queue[founds++] = dest;
							BitMask[dest] = 1;
						}
						else if (destination.y != current.y && destination.y < MAX_SIZE){	
							Matrix[ (current.y     * MAX_SIZE) + destination.y ] = true;	
							Matrix[ (destination.y * MAX_SIZE) + current.y 	   ] = true;	
						}
					}
				}
			}

			Write(SMemF1, &F2SizePtr[0], Queue, founds);
	
			level++;
			FrontierSize = F2SizePtr[0];

			GlobalWrite(counter, &GlobalCounter);

			__syncthreads();
			F2SizePtr[0] = 0;
		}
		//printf("Block %d has done his job!\n", Bid);
		FrontierSize = 1; 
		level = 0; 
		src = 0;
	}
}



/*
* KERNEL FOR BOTTOM-UP VISIT
*/
__global__ void BottomUp_Kernel(	const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									bool*  __restrict__ BitMask,
									const int BMsize,
									const int N)
{

	int founds = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH)
	{
		if(BitMask[index] == 0 && visit(devNode, devEdge, BitMask, index))
		{
			BitMask[index] = 1;
			founds++;
		}
	}
	GlobalWrite( founds, &BottomUp_FrontSize);
}