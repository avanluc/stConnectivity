#include <stConnectivity.hpp>

/*
* CUDA Kernel Device code
*/
__global__ void stConn(const int* devNodes, const int* devEdges, int2* Dist_Col, const int level, bool* newLevel, bool* Matrix, const int nof_nodes, const int nof_distNodes)
{

	// Calcolate thread id
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int id = tx + (bx*BLOCK_SIZE);
    //printf("Visita partita da %d\n",id);

	for (int i = id; i < nof_nodes; i += gridDim.x * BLOCK_SIZE) 
	{
	    if(i < nof_nodes)
	    {
	    	//printf("Visita partita da %d\n",i);
			if (Dist_Col[i].x == level) 										// se i è alla distanza giusta (ultima distanza visitata)
			{
				const int2 currDC = Dist_Col[i];
	    		//printf("Visita partita da %d con Distanza %d e Colore %d\n", i, currDC.x,  currDC.y);
	    		//printf("Visita partita da %d con Distanza %d e Colore %d\n", i,  Dist_Col[i].x,   Dist_Col[i].y);
	    		for (int j = devNodes[i]; j < devNodes[i + 1]; j++) 			// per ogni vicino di i
				{
					const int dest = devEdges[j];
					int2 destDC = Dist_Col[dest];
					//printf("dal nodo %d vedo il nodo %d e destDC.x != INT_MAX vale %d\n", i, dest, (destDC.x != INT_MAX));

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
						Dist_Col[dest].x = level + 1;							// lo visito e aggiorno la distanza
						Dist_Col[dest].y = currDC.y;						// setto il suo Dist_Cole come il Dist_Cole del padre
						newLevel[0] = true;										// notifico che ho trovato nuovi nodi da visitare
					}
				}
			}
		}
    }

}
