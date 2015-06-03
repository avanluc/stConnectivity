#include <stConnectivity.hpp>

#define	BLOCK_FRONTIER_LIMIT 1024
#define	     SMEMORY_SIZE	1024
#define		   BLOCK_SIZE	1024
#define		Thread_Per_SM	2048
#define          N_OF_SMs	12
#define    MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)
#define    MAX_CONCURR_BL(BlockDim)	( MAX_CONCURR_TH / (BlockDim) )
#define 			  Tid 	threadIdx.x


__device__ unsigned int GSync[MAX_CONCURR_BL(BLOCK_SIZE)];
__device__ bool devNextLevel[4];


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


__device__ __forceinline__ int LaneID() {
    int ret;
    asm("mov.u32 %0, %laneid;" : "=r"(ret) );
    return ret;
}


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


/*
*	CUDA Kernel core
*/
__device__  __forceinline__ bool stConnCore(const int* __restrict__ devNodes, 
											const int* __restrict__ devEdges, 
											int2* __restrict__ Dist_Col, 
											const int nof_nodes, 
											const int nof_distNodes,
											const int id, 
											const int level, 
											bool* Matrix){
	bool newLevel = false;
	for (int i = id; i < nof_nodes; i += gridDim.x * BLOCK_SIZE) 
	{
	    if(i < nof_nodes)
	    {
			if (Dist_Col[i].x == level) 									// se i è alla distanza giusta (ultima distanza visitata)
			{
	    		const int2 currDC = Dist_Col[i];
	    		for (int j = devNodes[i]; j < devNodes[i + 1]; j++) 		// per ogni vicino di i
				{
					const int dest = devEdges[j];
					const int2 destDC = Dist_Col[dest];
					
					if(destDC.x != INT_MAX && destDC.y != currDC.y)			// se è già stato visitato da qualcun altro
					{
						Matrix[ currDC.y*nof_distNodes + destDC.y ] = true;	// aggiorna la matrice di adiacenza
						Matrix[ destDC.y*nof_distNodes + currDC.y ] = true;	// aggiorna la matrice di adiacenza
					}
					if (destDC.x == INT_MAX) 								// se non è ancora stato visitato
					{
						int2 val = make_int2(level + 1,currDC.y);
						atomicStore(&Dist_Col[dest], val);
						// Dist_Col[dest].x = level + 1;
						// Dist_Col[dest].y = currDC.y;
						newLevel = true;									// notifico che ho trovato nuovi nodi da visitare
					}
				}
			}
		}
    }
    return newLevel;
}



/*
* CUDA Kernel device code
*/
__global__ void stConn(	const int* __restrict__ devNodes, 
						const int* __restrict__ devEdges, 
						int2* __restrict__ Dist_Col, 
						const int nof_nodes, 
						const int nof_distNodes,
						bool* Matrix){

	// Calcolate thread id
    const int id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    int level = 0;
    
    do{
    	// Set NextLevel to false for the next iteration
	    if (id == 0)
			devNextLevel[(level +1) & 3] = false;
	    // Call Kernel core function
	    int newLevel = stConnCore(devNodes, devEdges, Dist_Col, nof_nodes, nof_distNodes, id, level, Matrix);
		
		// Check if any thread of each warp had set newLevel to true
		// And in that case only one thread per warp writes in global memory
	    if (__any(newLevel) && LaneID() == 0)
			devNextLevel[level & 3] = true;

    	level++;

    	GlobalSync();
    // Iterate while there is no newLevel
    }while(devNextLevel[(level-1) & 3] );
}








/*
* ############### OLD CUDA Kernel Device code #####################
*/
/*
__global__ void stConn1(	const int* __restrict__ devNodes, 
						const int* __restrict__ devEdges, 
						int2* __restrict__ Dist_Col, 
						int level, 
						const int nof_nodes, 
						const int nof_distNodes,
						bool* Matrix){
//						bool* newLevel) { 
	// Calcolate thread id
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int id = tx + (bx*BLOCK_SIZE);
    //do{
	    if (id == 0)
			devNextLevel[((level & 1) + 1) & 1] = false;
		int newLevel = false;
	    
		for (int i = id; i < nof_nodes; i += gridDim.x * BLOCK_SIZE) 
		{
		    if(i < nof_nodes)
		    {
				if (Dist_Col[i].x == level) 										// se i è alla distanza giusta (ultima distanza visitata)
				{
		    		// printf("Visita partita da %d\n",i);
					const int2 currDC = Dist_Col[i];
		    		for (int j = devNodes[i]; j < devNodes[i + 1]; j++) 			// per ogni vicino di i
					{
						const int dest = devEdges[j];
						int2 destDC = Dist_Col[dest];
						
						if(destDC.x != INT_MAX && destDC.y != currDC.y)// se è già stato visitato da qualcun altro
						{
							// printf("incontro tra il nodo %d (visitato da %d) e %d (visitato da %d)\n", i, currDC.y, dest, destDC.y);
							Matrix[ currDC.y*nof_distNodes + destDC.y ] = true;// aggiorna la matrice di adiacenza
							Matrix[ destDC.y*nof_distNodes + currDC.y ] = true;// aggiorna la matrice di adiacenza
						}
						if (destDC.x == INT_MAX) 							// se non è ancora stato visitato
						{
							// printf("dal nodo %d visito il nodo %d\n", i, dest);
							int2 val = make_int2(level + 1,currDC.y);
							atomicStore(&Dist_Col[dest], val);
							newLevel = true;										// notifico che ho trovato nuovi nodi da visitare
						}
					}
				}
			}
	    }
	    if (__any(newLevel) && LaneID() == 0)		// (id & 31) === LaneID()
	    {
			devNextLevel[level & 1] = true;
		//	printf("at level %d devNextLevel[%d] = %d\n",level, (level & 1), devNextLevel[level & 1]);
	    }
    	//if(id == 0){
    	//level++;
		//	printf("at level %d devNextLevel[%d] = %d\n",level, ((level-1) & 1), devNextLevel[(level-1) & 1]);
    	//}
    	//GlobalSync(id);
    //}while(devNextLevel[(level-1) & 1]);
}*/




/*__device__  __forceinline__ void MVKernel_gm(const bool* A, const bool * X, bool* Y, const int N)
{
	int bx = blockIdx.x; 
	int tx = threadIdx.x; 
  	int Row = bx * BLOCK_SIZE + tx;
  
	float Pvalue = 0;

	for (unsigned int k = 0; k < N; k++)         
		Pvalue += A[Row*N + k] * X[k];

	__syncthreads();

	if(Row < N)  		
		Y[Row] = Pvalue;
	__syncthreads();
}


__global__ void MatrixBFS(const bool* A, const bool * X, const int N, bool * Y)
{
	int bx = blockIdx.x; 
	int tx = threadIdx.x; 
  	int i = bx * BLOCK_SIZE + tx;
  	bool Union = 0;

  	//while()
  	//{
	  	if(i < N)
	  	{
		  	MVKernel_gm(A, X, Y, N);
		  	Union = Union & X[i];
		  	//X[i] = Y[i] & (Union==false);
		  		
	  	}
	  	GlobalSync();
  	//}
}*/


__global__ void BFS_BlockKernel (	const int* __restrict__	devNode,
									const int* __restrict__	devEdge,
									int* __restrict__	devDistance,
									const int* __restrict__	devSource,
									const int Nsources,
									const int N) {

	int Queue[N];
	int level = 0;
	int FrontierSize = Nsources;

	__shared__ int SMemF1[SMEMORY_SIZE];
	__shared__ int SMemF2[SMEMORY_SIZE];

	SMemF1 = devSource;

	while (FrontierSize && FrontierSize < BLOCK_FRONTIER_LIMIT) {

		int founds = 0;
		for (int t = Tid; t < FrontierSize; t += BLOCK_SIZE) {
			const int index = SMemF1[t];
			const int start = devNode[index];
			int end = devNode[index + 1];

			for (int k = start; k < end; k++) {
				const int dest = devEdge[k];

				if (devDistance[dest] == -1) {
					devDistance[dest] = level;
					Queue[founds++] = dest;
				}
			}
		}

		//int WarpPos, n, total;
		//singleblockQueueAdd(founds, F2SizePtr, WarpPos, n, total, level, (int*) &SMem[TEMP_POS]); 	//  Util/GlobalWrite.cu

		//swapDev(SMemF1, SMemF2);
		level++;
		__syncthreads();
		//FrontierSize = F2SizePtr[0];
	}
}
