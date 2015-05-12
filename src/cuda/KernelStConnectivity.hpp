#include <stConnectivity.hpp>

/*
* CUDA Kernel Device code
*/
// __global__ void stConn(const int* devNodes, const int* devEdges, int2* Dist_Col, const int level, bool* newLevel, bool* Matrix, const int nof_nodes, const int nof_distNodes)
// {
__global__ void stConn(	const int* __restrict__ devNodes, 
						const int* __restrict__ devEdges, 
						int2* __restrict__ Dist_Col, 
						const int level, 
						const int nof_nodes, 
						const int nof_distNodes,
						bool* Matrix){
//						bool* newLevel) { 

	// Calcolate thread id
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int id = tx + (bx*BLOCK_SIZE);
    if (id == 0)
		devNextLevel[((level & 1) + 1) & 1] = false;
	int newLevel = false;
    
	for (int i = id; i < nof_nodes; i += gridDim.x * BLOCK_SIZE) 
	{
	    if(i < nof_nodes)
	    {
			if (Dist_Col[i].x == level) 										// se i è alla distanza giusta (ultima distanza visitata)
			{
	    		//printf("Visita partita da %d\n",i);
				const int2 currDC = Dist_Col[i];
	    		for (int j = devNodes[i]; j < devNodes[i + 1]; j++) 			// per ogni vicino di i
				{
					const int dest = devEdges[j];
					int2 destDC = Dist_Col[dest];
					
					if(destDC.x != INT_MAX && destDC.y != currDC.y)// se è già stato visitato da qualcun altro
					{
						//printf("incontro tra il nodo %d (visitato da %d) e %d (visitato da %d)\n", i, currDC.y, dest, destDC.y);
						//printf("aggiorno la matrice in posizione %d e %d\n", (currDC.y*nof_distNodes + destDC.y), (destDC.y*nof_distNodes + currDC.y));
						Matrix[ currDC.y*nof_distNodes + destDC.y ] = true;// aggiorna la matrice di adiacenza
						Matrix[ destDC.y*nof_distNodes + currDC.y ] = true;// aggiorna la matrice di adiacenza
					}
					if (destDC.x == INT_MAX) 							// se non è ancora stato visitato
					{
						//printf("dal nodo %d visito il nodo %d\n", i, dest);
						int2 val = make_int2(level + 1,currDC.y);
						atomicStore(&Dist_Col[dest], val);
						// Dist_Col[dest].x = level + 1;							// lo visito e aggiorno la distanza
						// Dist_Col[dest].y = currDC.y;							// setto il suo Dist_Cole come il Dist_Cole del padre
						newLevel = true;										// notifico che ho trovato nuovi nodi da visitare
					}
				}
			}
		}
    }
    if (__any(newLevel) && LaneID() == 0)		// (Tid & 31) === LaneID()
    {
		devNextLevel[level & 1] = true;
		//printf("devNextLevel[%d] = %d\n", (level & 1), devNextLevel[level & 1]);
    }
}

/*
__global__ void Mat(	const int* __restrict__ devNodes, 
						const int* __restrict__ devEdges, 
						int2* __restrict__ Dist_Col, 
						const int level, 
						const int nof_nodes, 
						const int nof_distNodes,
						bool* Matrix,
						bool* newLevel) { 
for (int i = ID; i < nof_nodes; i += gridDim.x * BlockDim) {
		if (Distance[i] == level) {
			for (int j = devNodes[ID]; j < devNodes[ID + 1]; j++) {
				const int dest = devEdges[j];
				if (Distance[dest] == INT_MAX) {
					Distance[dest] = level + 1;
					newLevel = true;
				}
			}
		}
	}
}
*/