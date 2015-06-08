
#include "prefixSumAsm.cu"

#define	BLOCK_FRONTIER_LIMIT 1024
#define	     SMEMORY_SIZE	1024
#define		   BLOCK_SIZE	1024
#define		Thread_Per_SM	2048
#define          N_OF_SMs	12
#define    MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)
#define    MAX_CONCURR_BL(BlockDim)	( MAX_CONCURR_TH / (BlockDim) )
#define 			  Tid 	threadIdx.x

/*#define HASHTABLE_BLOCK_POS  0
#define END_OF_HASHTABLE	(4096 * 8)	// 8: long long int size
#define       F2Size_POS	END_OF_HASHTABLE
#define         TEMP_POS	(F2Size_POS + 4)
#define     END_TEMP_POS	(TEMP_POS + 32 * 4)
#define    FRONTIER_SIZE	(((49152 - END_TEMP_POS) / 2) - 2)		//-2 align
#define     F1_BLOCK_POS	(END_TEMP_POS)
#define     F2_BLOCK_POS	(F1_BLOCK_POS + FRONTIER_SIZE)
*/

#define SM_BYTE_PER_BLOCK	40960
#define SM_BIT_PER_BLOCK	SM_BYTE_PER_BLOCK * 8
#define         F1_OFFSET	0
#define         F2_OFFSET	SM_BIT_PER_BLOCK/2


__device__ unsigned int GSync[MAX_CONCURR_BL(BLOCK_SIZE)];
__device__ bool devNextLevel[4];

__device__ int GlobalCounter;

extern __shared__ unsigned char SMem[];

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


/*__device__ __forceinline__ int LaneID() {
    int ret;
    asm("mov.u32 %0, %laneid;" : "=r"(ret) );
    return ret;
}*/


__global__ void GReset() {
	int bx = blockIdx.x;
    int tx = threadIdx.x;
    int id = tx + (bx*BLOCK_SIZE);
	if (id < MAX_CONCURR_BL(BLOCK_SIZE))
		GSync[id] = 0;
}


/*
*	Glocal Syncronization function
*/
__device__  __forceinline__ void GlobalSync() {
	volatile unsigned *VolatilePtr = GSync;
	__syncthreads();
	
	if (blockIdx.x == 0) {
		if (threadIdx.x == 0){
			VolatilePtr[blockIdx.x] = 1;
		}
		//__syncthreads();

		if (threadIdx.x < MAX_CONCURR_BL(BLOCK_SIZE))
			while ( cub::ThreadLoad<cub::LOAD_CG>(GSync + threadIdx.x) == 0 );

		__syncthreads();

		if (threadIdx.x < MAX_CONCURR_BL(BLOCK_SIZE)){
			VolatilePtr[threadIdx.x] = 0;
		}
	}
	else {
		if (threadIdx.x == 0) {
			VolatilePtr[blockIdx.x] = 1;
			while (cub::ThreadLoad<cub::LOAD_CG>(GSync + blockIdx.x) == 1);
		}
		__syncthreads();
	}
}


__device__ __forceinline__ void FrontierReserve_Warp(int* GlobalCounter, int founds, int& n, int &totalWarp, int& globalOffset) {
		n = founds;
		totalWarp = warpExclusiveScan<32>(n);
		int oldCounter;
		if (LaneID() == 0)	// && totale != 0)
			oldCounter = atomicAdd(GlobalCounter, totalWarp);

		globalOffset = __shfl(oldCounter, 0);
}



__device__ __forceinline__ void Write(int* devFrontier, int* GlobalCounter, int* Queue, int founds) {
		
		int n, total, globalOffset;
		FrontierReserve_Warp(GlobalCounter, founds, n, total, globalOffset);

		const int pos = globalOffset + n;
		for (int i = 0; i < founds; i++){
			devFrontier[pos + i] = Queue[i];
			printf("thread %d writing %d at position %d\n", Tid, Queue[i], pos + i);
		}
}

__device__ __forceinline__ void swapDev(int*& A, int*& B) {
	int* temp = A;	// frontiers swap
	A = B;
	B = temp;
}

__global__ void BFS_BlockKernel (	const int* __restrict__	devNode,
									const int* __restrict__	devEdge,
									int2* __restrict__	devDistance,
									int* __restrict__	devSource,
									const int Nsources) {
	int Queue[REG_QUEUE];
	int level = 0;
	int FrontierSize = Nsources;

	int* SMemF1 = (int*) &SMem[F1_OFFSET];
	int* SMemF2 = (int*) &SMem[F2_OFFSET];

	if (Tid < FrontierSize)
		SMemF1[Tid] = devSource[Tid]; 

	//while (FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT) {

		int founds = 0;
		for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE) {
			const int index = SMemF1[t];
			const int start = devNode[index];
			int end = devNode[index + 1];

			for (int k = start; k < end; k++) {
				const int dest = devEdge[k];

				if (devDistance[dest].x == INT_MAX) {
					int2 temp = make_int2(level, 0);
					atomicStore(&devDistance[dest], temp);
					//devDistance[dest].x = level;
					Queue[founds++] = dest;
				}
			}
		}

		int WarpPos, n, total;
		Write(SMemF2, &GlobalCounter, Queue, founds); 	//  Util/GlobalWrite.cu

		swapDev(SMemF1, SMemF2);
		level++;
		__syncthreads();
		FrontierSize = GlobalCounter;
	//}
}
