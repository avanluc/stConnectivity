#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"


__device__ int GlobalCounter = 0;
__device__ int exitFlag = 0;
__device__ int connected = 0;
__device__ int GQueue[REG_QUEUE];


/*
* KERNEL FUNCTION FOR STCONNECTIVITY
*/
__global__ void STCONN_BlockKernel (const int* __restrict__ devNode,
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



/*
* KERNEL FUNCTION FOR BFS
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



/*
* PARALLEL MATRIX BFS !!!DEPRECATED!!!
*/
__global__ void MatrixBFS1(	bool* __restrict__ Matrix,
							int* __restrict__ Visited,
							const int src,
							const int dest,
							const int Nsources)
{
	int  Queue[REG_QUEUE];
	int  FrontierSize = 1;
	int* SMemF1 	= (int*) &SMem[F1_OFFSET];
	int* F2SizePtr 	= (int*) &SMem[F2Size_POS];
	int* Block_Exit = (int*) &SMem[EXIT_FLAG];
	int* SMemTemp 	= (int*) &SMem[TEMP_POS];
	int founds = 0;


	if(GTid < Nsources)
	{

		if(GTid == 0){
			F2SizePtr[0] = 0;
			Visited[src] = 1;		
		}

		/***    PRIMA FASE: VISITA SRC    ***/
		if( Matrix[ src * Nsources + GTid ] )
		{
			//printf("GTID = %d\n", GTid);
			Visited[GTid] = 1;
			Queue[founds++] = GTid;
		}

		Write(SMemF1, &F2SizePtr[0], Queue, Block_Exit, founds);

		__syncthreads();
		FrontierSize = F2SizePtr[0];
		//if(Tid == 0)
			//printf("FrontierSize%d = %d, SMEM = %d, %d\n",Bid, FrontierSize, SMemF1[0], SMemF1[1] );
		__syncthreads();
		F2SizePtr[0] = 0;


		/***    SECONDA FASE: VISITA MATRICE    ***/
		while(FrontierSize)
		{
			founds = 0;

			if(Tid < 34){
				SMemTemp[Tid] = 0;
				//printf("Warp %d = %d\n", Tid, SMemTemp[Tid]);
			}


			for (int i = 0; i < FrontierSize; ++i)
			{
				int qNode = SMemF1[i];
					//if(Tid == 0)
					//	printf("start visit from %d = %d\n", qNode, SMemF1[i] );
					for (int j = Tid; j < Nsources; j+=BLOCK_SIZE)
					{
						if( Matrix[ qNode * Nsources + j ] && Visited[j] == 0)
						{
							//printf("\tThread %d aggiunge %d\n", Tid, j );
							Visited[j] = 1;
							Queue[founds++] = j;
						}
					}
			}

			Write(SMemF1, &F2SizePtr[0], Queue, Block_Exit, founds);
			__syncthreads();
			FrontierSize = F2SizePtr[0];
			//if(Tid == 0)
			//	printf("FrontierSize%d = %d\n",Bid, FrontierSize );
			__syncthreads();
			F2SizePtr[0] = 0;
		}
		
		//__syncthreads();

		if(Tid == 0)
			connected = Visited[dest];
	}
}


