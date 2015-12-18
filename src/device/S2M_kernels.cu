#include <assert.h>
#include <curand_kernel.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "S2M_kernel_cores.cu"

__device__ int BottomUp_FrontSize  = 0;
__device__ int BottomUp_FrontSize1 = 0;
__device__ int FrontierOffset = 1;
__device__ int GlobalCounter = 0;
__device__ int VisitResult  = 0;
__device__ int VisitResult1 = 0;
__device__ int color = 0;
__device__ bool matrixResult = 0;


/*
* INITIALIZATION OF CURAND SEED
*/
__global__ void setup_curand(curandState* __restrict__ state)
{
	if(Tid == 0)
		curand_init(Bid + clock(), 0, 0, &state[Bid]);
}



/*
* INITIALIZATION OF CURAND SEED
*/
__global__ void init_bottomUp_only(const int* __restrict__ devSource, bool* __restrict__ BitMask)
{
	if(GTid < 2)
		BitMask[devSource[GTid]] = 1;
}



/*
* RESET OF SHARED MEMORY
*/
template <int BlockDim>
__global__ void clean()
{
	if(Tid < IntSMem_Per_Block(BlockDim))
		SMem[Tid] = 0;
}



/*
* CHECK IF EVERY VERTEX IS VISITED
*/
__global__ void CheckVisit(	const bool* __restrict__ BitMask, const int2* __restrict__ devDistance, const int N)
{
	int founds = 0;
	int founds1 = 0;
	for (int index = GTid; index < N; index += MAX_CONCURR_TH)
	{
		if(BitMask[index] == 0)
		{
			founds++;
		}
		if(devDistance[index].y == INT_MAX)
		{
			founds1++;
		}
	}
	atomicAdd(&VisitResult, founds);
	atomicAdd(&VisitResult1, founds1);
}



/*
* KERNEL FUNCTION FOR STCONNECTIVITY
*/
template <int BlockDim, int WARP_SZ>
__global__ void TopDown_Kernel (const int* __restrict__ devNode,
								const int* __restrict__ devEdge,
								const int* __restrict__ devSource,
								int2* __restrict__ devDistance,
								bool* __restrict__ Matrix,
								bool* __restrict__ BitMask, 
								curandState*  __restrict__ globalState,
								const double Treshold,
								const int E,
								const int N)
{
	int Queue[REG_QUEUE];
	int FrontierSize = 1;
	int level = 1;
	uint src  = 0;
	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* F2SizePtr = (int*) &SMem[F2Size_POS];

	const int VirtualID = Tid / WARP_SZ;
	const int Stride =  BlockDim / WARP_SZ;


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
			for (int t = VirtualID; t < FrontierSize; t += Stride)
			{
				const int index = SMemF1[t];
				const int start = devNode[index];
				const int   end = devNode[index + 1];
				const int2 current = devDistance[index];
	
				for (int k = start + (Tid % WARP_SZ); k < end; k+=WARP_SZ)
				{
					edgeVisit(devEdge[k], current, devDistance[devEdge[k]], level, Queue, devDistance, BitMask, Matrix, founds, counter);
				}
			}

			Write(SMemF1, &F2SizePtr[0], Queue, founds);
			
			level++;
			FrontierSize = F2SizePtr[0];

			AtomicWrite(counter, &GlobalCounter);

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
	AtomicWrite( founds, &BottomUp_FrontSize );
}




/*
* KERNEL FOR ADJACENCE MATRIX VISIT
*/
__global__ void matrixVisit(const bool* __restrict__ Matrix,
							bool* __restrict__ BitMask,
							int* __restrict__ Queue,
							const int totSources,
							const int src,
							const int dst)
{
	int lQueue[REG_QUEUE];
	int oldFrontierSize = 1;
	int founds = 0;
	int qPos = 0;
	Queue[qPos] = src;
	int level = 0;

	if(GTid < totSources)
	{
		while(level < 2 /*FrontierOffset || !BitMask[dst]*/)
		{
			if(GTid == 0)
				printf("Visiting row %d\n", Queue[qPos]);

			if(Matrix[ (Queue[qPos]*MAX_SIZE) + GTid ] && !BitMask[GTid])
			{
				BitMask[GTid] = 1;
				lQueue[founds++] = GTid;
			}

			Write(&Queue[oldFrontierSize], &FrontierOffset, lQueue, founds);

			__syncthreads();
			oldFrontierSize += FrontierOffset;
			founds = 0;
			qPos++;
			level++;
		}
	}

	if(GTid == 0)
		matrixResult = BitMask[dst];
}