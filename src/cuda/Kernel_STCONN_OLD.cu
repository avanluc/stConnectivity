#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "CacheFunc.cu"



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
__global__ void STCONN_BlockKernel (const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									const int* __restrict__ devSource,
									int2* __restrict__ devDistance,
									bool* __restrict__ Matrix,
									bool*  __restrict__ BitMask, 
									const int Nsources,
									const int E)
{
	int Queue[REG_QUEUE];
	int FrontierSize = 1;
	int level = 0;
	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];

	for(int j = Bid; j < Nsources; j += gridDim.x)
	{
		if (Tid < FrontierSize){
			SMemF1[Tid] = devSource[j];
			BitMask[SMemF1[Tid]] = 1;
		}
	
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
				}
			}

			Write(SMemF1, &F2SizePtr[0], Queue, founds);
	
			level++;
			FrontierSize = F2SizePtr[0];

			GlobalWrite(counter, &GlobalCounter);

			__syncthreads();
			F2SizePtr[0] = 0;

			if(__int2double_rn(GlobalCounter) / __int2double_rn(E) > TRESHOLD)
				return;
		}
		level = 0; FrontierSize = 1;
	}
}



/*
* KERNEL FOR BOTTOM-UP VISIT
*/
__global__ void Bottom_Up_Kernel(	const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									bool*  __restrict__ BitMask,
									const int BMsize,
									const int N)
{

	/*const int stride = gridDim.x * BLOCK_SIZE * 8;
	bool *BitMarkArray = BitMask + Tid * 4;
	int founds  = 0;
	for (int BlockIndex = Bid * BLOCK_SIZE * 8; BlockIndex < N; BlockIndex += stride)
	{
		bool Queue[8];
 
		reinterpret_cast<int*>(Queue)[0] = reinterpret_cast<int*>(BitMask + BlockIndex)[Tid];
		reinterpret_cast<int*>(Queue)[1] = __ldg( &reinterpret_cast<int*>(BitMask + BlockIndex)[Tid + BLOCK_SIZE] );
		
		Queue[0] = BitMask[BlockIndex + (Tid * 4) + 0];
		Queue[1] = BitMask[BlockIndex + (Tid * 4) + 1];
		Queue[2] = BitMask[BlockIndex + (Tid * 4) + 2];
		Queue[3] = BitMask[BlockIndex + (Tid * 4) + 3];
		Queue[4] = BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 0];
		Queue[5] = BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 1];
		Queue[6] = BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 2];
		Queue[7] = BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 3];

		#pragma unroll
		for (int i = 0; i < 8; i++){
			const int ldg_stride = i >= 4 ? BLOCK_SIZE * 4 : 0;
			const int index = BlockIndex + (Tid * 4) + ldg_stride + i%4;
			if (Queue[i] == 0 && index < N && visitAdjiacent(index, devNode, devEdge, BitMask))
			{
				Queue[i] = 1;
				founds++;
			}
		}

		BitMask[BlockIndex + (Tid * 4) + 0] = Queue[0];
		BitMask[BlockIndex + (Tid * 4) + 1] = Queue[1];
		BitMask[BlockIndex + (Tid * 4) + 2] = Queue[2];
		BitMask[BlockIndex + (Tid * 4) + 3] = Queue[3];
		BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 0] = Queue[4];
		BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 1] = Queue[5];
		BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 2] = Queue[6];
		BitMask[BlockIndex + (Tid * 4) + (BLOCK_SIZE * 4) + 3] = Queue[7];

		reinterpret_cast<int*>(BitMask + BlockIndex)[Tid] = reinterpret_cast<int*>(Queue)[0];
		reinterpret_cast<int*>(BitMask + BlockIndex)[Tid + BLOCK_SIZE] = reinterpret_cast<int*>(Queue)[1];
		
		BitMarkArray += stride;
	}
	GlobalWrite( founds, &BottomUp_FrontSize);*/




	int founds = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH)
		if(BitMask[index] == 0 && visitAdjiacent(devNode, devEdge, BitMask, index))
		{
			BitMask[index] = 1;
			founds++;
		}
	GlobalWrite( founds, &BottomUp_FrontSize);
}