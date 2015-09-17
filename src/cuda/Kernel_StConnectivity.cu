#pragma once

#include <assert.h>
#include "prefixSumAsm.cu"
#include "GlobalSync.cu"
#include "definition.cuh"

__device__ int GlobalCounter = 0;
__device__ int exitFlag = 0;
__device__ int connected = 0;
__device__ int GQueue[REG_QUEUE];

extern __shared__ unsigned char SMem[];


 /*
* Assert for CUDA functions
*/
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}



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


/*
* 
*/
__device__ __forceinline__ void FrontierReserve_Block(int* Front_size, int founds, int& n, int &totalBlock, int& globalOffset){
	int* SM = (int*) &SMem[TEMP_POS];
	n = founds;
	const int warpId = WarpID();
	SM[warpId] = warpExclusiveScan<32>(n);
	__syncthreads();
	if (Tid < BLOCK_SIZE / 32)
	{
		int sum = SM[Tid];
		const int total = warpExclusiveScan<BLOCK_SIZE / 32>(sum);

		if (Tid == 0)
		{
			SM[32] = total;
			SM[33] = atomicAdd(Front_size, total);
		}
		SM[Tid] = sum;
	}
	__syncthreads();

	n += SM[warpId];
	totalBlock = SM[32];
	globalOffset = SM[33];
}


/*
*
*/
__device__ __forceinline__ void Write(int* devFrontier, int* Front_size, int* Queue, int* Exit_Flag, int founds) {
		
	int n, total, globalOffset;
	FrontierReserve_Block(Front_size, founds, n, total, globalOffset);

	const int pos = globalOffset + n;
	for (int i = 0; i < founds; i++)
	{
		if((pos + i) < BLOCK_FRONTIER_LIMIT)
		{
			devFrontier[pos + i] = Queue[i];
		}
		else
		{
			#if(ATOMIC)
				Exit_Flag[0] = 1;
			#else
				exitFlag = 1;
			#endif
		}
	}
}


/*
*
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
	
		while (FrontierSize)
		{
			if(!ATOMIC && exitFlag)
				break;
			if(Block_Exit[0]){
				break;
			}
	
			int founds = 0; //int counter = 0;
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
	
					#if (!BFS && ATOMIC)
						
						//counter++;
						if(founds < REG_QUEUE)
						{
							if ( atomicCAS(&devDistance[dest].x, INT_MAX, level) == INT_MAX ) {
								devDistance[dest].y = current.y;
								Queue[founds++] = dest;
							}
							else if (destination.y != current.y && destination.y < Nsources){	
								Matrix[ (current.y     * Nsources) + destination.y ] = true;	
								Matrix[ (destination.y * Nsources) + current.y 	   ] = true;	
							}
						}
					
					#elif (BFS && ATOMIC)	
						//counter++;
						if(founds < REG_QUEUE)
						{
							if ( atomicCAS(&devDistance[dest].x, INT_MAX, level) == INT_MAX ) {
								devDistance[dest].y = current.y;
								Queue[founds++] = dest;
							}
						}
	
					#elif (BFS && !ATOMIC)
						if (destination.x == INT_MAX){	
							devDistance[dest].x = level;
							devDistance[dest].y = current.y;
							Queue[founds++] = dest;
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
				//atomicAdd(&GlobalCounter, counter);
				if (Tid == 0)
					printf("Visited %d\n",atomicAdd(&GlobalCounter, FrontierSize));
			#endif
	
			__syncthreads();
			F2SizePtr[0] = 0;
		}
		level = 0;
		FrontierSize = 1;
	}
}



__global__ void MatrixBFS1(	bool* __restrict__ Matrix,
							int* __restrict__ Visited,
							const int src,
							const int dest,
							const int Nsources)
{
	if(GTid < Nsources)
	{
		int  Queue[REG_QUEUE];
		int  FrontierSize = 1;
		int* SMemF1 	= (int*) &SMem[F1_OFFSET];
		int* F2SizePtr 	= (int*) &SMem[F2Size_POS];
		int* Block_Exit = (int*) &SMem[EXIT_FLAG];
		int founds = 0;

		if(GTid == 0)
			Visited[src] = 1;
	
		/***    PRIMA FASE: VISITA SRC    ***/
		printf("GTid = %d, mat = %d\n",GTid , Matrix[ src * Nsources + GTid ]);
		if( Matrix[ src * Nsources + GTid ] )
		{
			printf("GTid = %d\n",GTid );
			Visited[GTid] = 1;
			Queue[founds++] = GTid;
		}

		Write(SMemF1, &F2SizePtr[0], Queue, Block_Exit, founds);

		__syncthreads();
		FrontierSize = F2SizePtr[0];
		__syncthreads();
		F2SizePtr[0] = 0;


		/***    SECONDA FASE: VISITA MATRICE    ***/
		while(FrontierSize)
		{
			founds = 0;
			for (int i = 0; i < FrontierSize; ++i)
			{
				for (int j = Tid; j < Nsources; j+=BLOCK_SIZE)
				{
					if( j < Nsources && Matrix[ i * Nsources + j ] && Visited[j] == 0)
					{
						Visited[j] = 1;
						Queue[founds++] = j;
					}
				}
			}

			Write(SMemF1, &F2SizePtr[0], Queue, Block_Exit, founds);
			__syncthreads();
			FrontierSize = F2SizePtr[0];
			__syncthreads();
			F2SizePtr[0] = 0;
		}
		
		//__syncthreads();

		if(Tid == 0)
			connected = Visited[dest];
	}
}


