#pragma once

#include <stConnectivity.hpp>
#include <../../cub/cub.cuh>

#define		Thread_Per_SM	2048
#define          N_OF_SMs	12
#define    MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)

#define    MAX_CONCURR_BL(BlockDim)	( MAX_CONCURR_TH / (BlockDim) )

__device__ unsigned int GSync[MAX_CONCURR_BL(BLOCK_SIZE)];

__global__ void GReset() {
	int bx = blockIdx.x;
    int tx = threadIdx.x;
    int id = tx + (bx*BLOCK_SIZE);
	if (id < MAX_CONCURR_BL(BLOCK_SIZE))
		GSync[id] = 0;
}

__device__  __forceinline__ void GlobalSync(int id) {
	volatile unsigned *VolatilePtr = GSync;
	__syncthreads();
	
	if (blockIdx.x == 0) {
		if (id == 0)
			VolatilePtr[blockIdx.x] = 1;
		//__syncthreads();

		if (id < MAX_CONCURR_BL(BLOCK_SIZE))
			while ( cub::ThreadLoad<cub::LOAD_CG>(GSync + id) == 0 );

		__syncthreads();

		if (id < MAX_CONCURR_BL(BLOCK_SIZE))
			VolatilePtr[id] = 0;
	}
	else {
		if (id == 0) {
			VolatilePtr[blockIdx.x] = 1;
			while (cub::ThreadLoad<cub::LOAD_CG>(GSync + blockIdx.x) == 1);
		}
		__syncthreads();
	}
}



__device__ __forceinline__ void memcpy_SIMD(int W_OFF, int cnt, int* dst, const int* src){
	for (int IDX = W_OFF; IDX < cnt; IDX += WARP_SIZE){
		dst[IDX] = src[IDX];
		printf("Copyng  nodes[%d] = %d = %d\n",IDX, src[IDX], dst[IDX]);
	}
	__threadfence_block();
}

__device__ __forceinline__ void memcpy_SIMD(int W_OFF, int cnt, int2* dst, int2* src){
	for (int IDX = W_OFF; IDX < cnt; IDX += WARP_SIZE){
		dst[IDX].x = src[IDX].x;
		dst[IDX].y = src[IDX].y;
		//printf("Copyng  Dist_Col[%d].x= %d = %d\n",IDX, src[IDX].x, dst[IDX].x);
	}
	__threadfence_block();
}

/*
__device__ __forceinline__ void expand_BFS_SIMD(int W_OFF, int cnt, const int* edges, int2* Dist_Col, const int level, bool newLevel){
	for (int IDX = W_OFF; IDX < cnt; IDX += WARP_SIZE){
		int v = edges[IDX];
		if(Dist_Col[v].x == INT_MAX){
			int2 val = make_int2(level + 1,currDC.y);
			atomicStore(&Dist_Col[dest], val);
			newLevel = true;
		}
		if(Dist_Col[v].x != INT_MAX && Dist_Col[v].y != currDC.y)// se è già stato visitato da qualcun altro
		{
			Matrix[ currDC.y*nof_distNodes + destDC.y ] = true;// aggiorna la matrice di adiacenza
			Matrix[ destDC.y*nof_distNodes + currDC.y ] = true;// aggiorna la matrice di adiacenza
		}
	}

	__threadfence_block();
}
*/

/*
* CUDA Kernel Device code
*/
// __global__ void stConn(const int* devNodes, const int* devEdges, int2* Dist_Col, const int level, bool* newLevel, bool* Matrix, const int nof_nodes, const int nof_distNodes)
// {
__global__ void stConn(	const int* __restrict__ devNodes, 
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
  //   int W_OFF = id % WARP_SIZE;
  //   int W_ID = id / WARP_SIZE;
   
  //   extern __shared__ warpmem_t SMEM[];
  //   warpmem_t *MY = SMEM + (tx / WARP_SIZE);

  //   if(id < nof_nodes){

	 //    int start = W_ID*BLOCK_SIZE;
	 //    int cnt = min(BLOCK_SIZE, nof_nodes-start);
	 //    // launch memcpy_SIMD with counter = min(BLOCK_SIZE, nof_nodes-start)
	 //    //memcpy_SIMD(W_OFF, cnt + 1, MY->nodes, &devNodes[start]);
	 //    for (int IDX = W_OFF; IDX < cnt+1; IDX += WARP_SIZE){
	 //    	MY->nodes[IDX] = devNodes[start+IDX];
		// 	printf("MY->nodes[%d] = %d\n", IDX, MY->nodes[IDX]);
	 //    }
		// __threadfence_block();
		// printf("MY->nodes[%d] = %d\n", W_OFF, MY->nodes[W_OFF]);
	 //    memcpy_SIMD(W_OFF, cnt, MY->Dst_Col, &Dist_Col[start]);

	 //    if (id == 0)
		// 	devNextLevel[((level & 1) + 1) & 1] = false;
		// int newLevel = false;

	 //    for (int v = 0; v < cnt; v++)
	 //    {
	 //    	if(MY->Dst_Col[v].x == level)
	 //    	{
	 //    		int num_vicini = MY->nodes[v+1] - MY->nodes[v];
	 //    		const int* vicini = &devEdges[MY->nodes[v]];
	 //    		//printf("MY->nodes[%d] - MY->nodes[%d] = %d\n", v+1, v, MY->nodes[v+1] - MY->nodes[v]);
	 //    		//printf("Visita partita da %d con grado %d\n", v, num_vicini);
	    		
	 //    		for (int IDX = W_OFF; IDX < num_vicini; IDX += WARP_SIZE)
	 //    		{
		// 			int n = vicini[IDX];
		// 			int2 nDC = Dist_Col[n];
		// 			if(nDC.x == INT_MAX){
		// 				int2 val = make_int2(level + 1,nDC.y);
		// 				atomicStore(&Dist_Col[n], val);
		// 				newLevel = true;
		// 			}
		// 			if(nDC.x != INT_MAX && nDC.y != MY->Dst_Col[v].y)// se è già stato visitato da qualcun altro
		// 			{
		// 				Matrix[ MY->Dst_Col[v].y*nof_distNodes + nDC.y ] = true;// aggiorna la matrice di adiacenza
		// 				Matrix[ nDC.y*nof_distNodes + MY->Dst_Col[v].y ] = true;// aggiorna la matrice di adiacenza
		// 			}
		// 		}
		// 		__threadfence_block();
		// 	}
	 //    }

	 //    if (__any(newLevel) && LaneID() == 0)		// (id & 31) === LaneID()
	 //    {
		// 	devNextLevel[level & 1] = true;
		// 	printf("at level %d devNextLevel[%d] = %d\n",level, (level & 1), devNextLevel[level & 1]);
	 //    }
		// memcpy_SIMD(W_OFF, min(BLOCK_SIZE, nof_nodes-start), &Dist_Col[start], MY->Dst_Col);    
  //   }
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
    	//__syncthreads();
    //}while(devNextLevel[(level-1) & 1]);
}





__global__ void Mat( bool* Matrix, int* Distance, int level, const int nof_nodes) { 

	const int ID = blockIdx.x * BLOCK_SIZE + threadIdx.x;

	do{
		if (ID == 0)
			devNextLevel[(level & 1) + 1] = false;
		int newLevel = false;

		for (int i = ID; i < nof_nodes; i += gridDim.x * BLOCK_SIZE) {
			if (Distance[i] == level) {
				//printf("Visita partita da %d\n",i);
				for (int j = 0; j < nof_nodes; j++) {
					//printf("Matrix[%d] = %d e Distance[%d] = %d\n",i*nof_nodes+j, Matrix[i*nof_nodes+j],j,Distance[j]);
					if (Matrix[i*nof_nodes + j] && Distance[j] == INT_MAX) {
						//printf("dal nodo %d visito il nodo %d\n", i, j);
						Distance[j] = level + 1;
						newLevel = true;
					}
				}
			}
		}
		if (__any(newLevel) && LaneID() == 0)		// (Tid & 31) === LaneID()
			devNextLevel[level & 1] = true;

		GlobalSync(threadIdx.x);

		if(ID == 0)
    		level++;

	}while(devNextLevel[(level-1) & 1] == true);
}