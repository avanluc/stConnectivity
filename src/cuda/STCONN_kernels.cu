#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "CacheFunc.cu"
#include "STCONN_kernel_cores.cu"
#include <curand_kernel.h>


/*
* INITIALIZATION OF CURAND SEED
*/
__global__ void setup_curand(curandState* state)
{
	if(Tid == 0)
		curand_init(Bid + clock(), 0, 0, &state[Bid]);
}



/*
* RESET OF SHARED MEMORY
*/
__global__ void clean()
{
	if(Tid < IntSMem_Per_Block(BLOCK_SIZE))
		SMem[Tid] = 0;
}



/*
* CHECK IF EVERY VERTEX IS VISITED
*/
__global__ void CheckVisit(	const bool* __restrict__ BitMask, const int2* __restrict__ devDistance, const int N)
{
	int founds = 0;
	int founds1 = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH){
		if(BitMask[index] == 0){
			founds++;
		}
		if(devDistance[index].y == INT_MAX){
			founds1++;
		}
	}
	atomicAdd(&VisitResult, founds);
	atomicAdd(&VisitResult1, founds1);
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
	int level = 1;
	uint src  = 0;
	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];

	if(Bid < 2 && Tid == 0)
	{
		src = devSource[Bid];
		BitMask[src] = 1;
		SMemF1[Tid] = src;
	}

	while(__int2double_rn(GlobalCounter) / __int2double_rn(E) < Treshold && color < MAX_SIZE )
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

		while ( FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT)
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
							BitMask[dest] = 1;
							Queue[founds++] = dest;
						}
						else if (destination.y < MAX_SIZE && destination.y != current.y ){
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
			if(__int2double_rn(GlobalCounter) / __int2double_rn(E) >= Treshold)
				return;
			
			F2SizePtr[0] = 0;
		}
		FrontierSize = 1; 
		level = 0; 
		src = 0;
	}
}



/*
* KERNEL FOR BOTTOM-UP VISIT
*/
__global__ void BottomUp_Kernel(const int* __restrict__ devNode,
								const int* __restrict__ devEdge,
								int2* __restrict__ devDistance,
								bool* __restrict__ Matrix,
								bool*  __restrict__ BitMask,
								const int N)
{
	int founds = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH)
	{
		if( BitMask[index] == 0 )
		{
			const int start = devNode[index];
			const int end = devNode[index + 1];
			for (int k = start; k < end; k++)
			{
				findConnection(devEdge[k], index, devDistance[index], devDistance[devEdge[k]], devDistance, BitMask, Matrix, founds);
			}
		}
	}
	GlobalWrite( founds, &BottomUp_FrontSize );
}










/*__global__ void BottomUp_Kernel(const int* __restrict__ devNode,
								const int* __restrict__ devEdge,
								int2* __restrict__ devDistance,
								bool* __restrict__ Matrix,
								bool*  __restrict__ BitMask,
								const int BMsize,
								const int N)
{
	int founds = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH)
	{
		if(BitMask[index] == 0)
		{
			const int start = devNode[index];
			const int end = devNode[index + 1];
			for (int k = start; k < end; k++)
			{
				const int dest = devEdge[k];
				const int2 destination = devDistance[dest];
				const int2 current = devDistance[index];
				if(BitMask[dest] == 1 )
				{
					if( BitMask[index] == 0 && destination.y != INT_MAX)
					{
						devDistance[index].x = destination.x + 1;
						devDistance[index].y = destination.y;
						BitMask[index] = 1;
						founds++;
					}
					else if( current.y != destination.y && current.y < MAX_SIZE && destination.y < MAX_SIZE )
					{
						Matrix[ (current.y 	   * MAX_SIZE) + destination.y ] = true;
						Matrix[ (destination.y * MAX_SIZE) + current.y     ] = true;
					}
				}
			}
		}
	}
	GlobalWrite( founds, &BottomUp_FrontSize);
}*/