#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "CacheFunc.cu"



/*
* kERNEL FUNCTION THAT CHECKS IF EVERY NODE HAS BEEN VISITED
*/
__global__ void CheckVisit(	const int* __restrict__ devNode, 
							const int* __restrict__ devEdge,
							int* __restrict__ BitMask,
							const int BMsize,
							const int N)
{
	int founds = 0;
	for (int i = GTid; i < N; i+=MAX_CONCURR_TH)
		if(markAccess<cub::LOAD_CS, int>(BitMask, i) == 0 ){
			founds++;
	}

	atomicAdd(&VisitResult, founds);
}


/*
* KERNEL FUNCTION FOR SET EXTRA BITS TO 1
*/
__global__ void initBitMask(int* BitMask, const int N, const int MaskSize)
{
	int totBits = MaskSize * 8 * sizeof(int);
	if(GTid < (totBits - N))
		markWrite<cub::LOAD_CS, cub::STORE_CS, int>(BitMask, N + GTid);
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
			markWrite<cub::LOAD_CG, cub::STORE_CG, int>(BitMask, SMemF1[Tid]);
			//AtomicMarkWrite(BitMask, SMemF1[Tid]);
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
							//AtomicMarkWrite(BitMask, dest);
							markWrite<cub::LOAD_CG, cub::STORE_CG, int>(BitMask, dest);
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
		FrontierSize = 1;  level = 0; 
	}
}



/*
* BOTTOM UP KERNEL
*/
__global__ void Bottom_Up_Kernel(	const int* __restrict__ devNode,
									const int* __restrict__ devEdge,
									int*  __restrict__ BitMask,
									const int BMsize,
									const int N)
{
	const int stride = gridDim.x * BLOCK_SIZE * 8;
	int* BitMarkArray = BitMask + Tid * 4;
	int founds  = 0;
	for (int BlockIndex = Bid * BLOCK_SIZE * 8; BlockIndex < BMsize; BlockIndex += stride)
	{
		char Queue[32];

		reinterpret_cast<int4*>(Queue)[0] = reinterpret_cast<int4*>(BitMarkArray + BlockIndex)[0];
		reinterpret_cast<int4*>(Queue)[1] = __ldg( &reinterpret_cast<int4*>(BitMarkArray + BlockIndex)[BLOCK_SIZE] );

		#pragma unroll
		for (int i = 0; i < 32; i++){
			#pragma unroll
			for (int j = 0; j < 8; j++){
				const int ldg_stride = i >= 16 ? BLOCK_SIZE * 128 : 0;
				const int index = (BlockIndex * 32) + (Tid * 128) + ldg_stride + (i%16 * 8) + j;
				if ((Queue[i] & (1 << (7 - j))) == 0 && index < N && visitAdjiacentBit(index, devNode, devEdge, BitMask))
				{
					Queue[i] |= (1 << (7 - j));
					founds++;
				}
			}
		}

		reinterpret_cast<int4*>(BitMarkArray + BlockIndex)[0] = reinterpret_cast<int4*>(Queue)[0];
		reinterpret_cast<int4*>(BitMarkArray + BlockIndex)[BLOCK_SIZE] = reinterpret_cast<int4*>(Queue)[1];
		
		BitMarkArray += stride;
	}
	GlobalWrite( founds, &BottomUp_FrontSize);
}
