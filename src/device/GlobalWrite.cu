#include "prefixSumAsm.cu"
#include "definition.cuh"

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
__device__ __forceinline__ void Write(int* devFrontier, int* Front_size, int* Queue, int founds) {
		
	int n, total, globalOffset;
	FrontierReserve_Block(Front_size, founds, n, total, globalOffset);

	const int pos = globalOffset + n;
	for (int i = 0; i < founds; i++)
	{
		if((pos + i) < BLOCK_FRONTIER_LIMIT)
		{
			devFrontier[pos + i] = Queue[i];
		}
	}
}



/*
* 
*/
__device__ __forceinline__ void FrontierReserve(int founds, int& n, int &totalBlock){
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
			SM[32] = total;
	}
	__syncthreads();
	totalBlock = SM[32];
}



/*
*
*/
__device__ __forceinline__ void AtomicWrite(int founds, int* GlobalVar) {
		
	int n, total;
	FrontierReserve(founds, n, total);
	
	if(Tid==0)
		atomicAdd(GlobalVar, total);
}
