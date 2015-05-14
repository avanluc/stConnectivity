#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <utility>
#include <climits>
#include <vector>
#include <math.h>
#include <bitset>

using namespace std;

// Block size
#define BLOCK_SIZE	256
#define WARP_SIZE	16
#define DIMGRID_MAX	65535

__device__ bool devNextLevel[2];

// Graph parameters (#nodes, #edges)
int  N, E;


// Edge structure
struct edge{
	int x;
	int y;
};

// Stuct for ordered nodes
struct ONodes{
	int id;
	int degree;
};

// Stuct for ordered nodes
struct warpmem_t{
	int nodes[BLOCK_SIZE+1];
	int2 Dst_Col[BLOCK_SIZE];
};


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
* Compare function for sorting
*/
int compare(const void *x, const void *y){
	edge a = *(edge *)x;
	edge b = *(edge *)y;
	return a.x < b.x ? -1 : (a.x > b.x ? 1 : (a.y < b.y ? -1 : a.y > b.y));
}

/*
* Compare function for sorting nodes per degree
*/
int ONodesCompare(const ONodes &a, const ONodes &b){
	return a.degree > b.degree;
}




bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
	vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato
	
	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata

	//printf("inizio bfs su matrice\n");

	while (left < right) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = 0; i < nof_nodes; ++i) 		// per ogni nodo nella lista di adiacenza di qNode
			if(adjMatrix[qNode*nof_nodes + i])
				if (!Visited[i]) {					// se dest non è ancora stato visitato 
					//printf("da %d visitato nuovo nodo %d\n", qNode, i);
					Visited[i] = true;				// lo marco come visitato
					Queue[right++] = i;				// lo aggiungo alla coda dei nodi incontrati durante la visita 
				}
	}
	delete [] Queue;
	return (Visited[target] == true);
}



vector< int > ChooseNodes(const int* Nodes, const int nof_nodes, const int nof_distNodes, const int source, const int target) {

	vector< int > distNodes(nof_distNodes);			// vector of distinguished nodes
	vector< ONodes > OrderedNodes(nof_nodes);		// vector of nodes ordered by degree

	distNodes[0] = source;
	
	if(nof_distNodes > 1)
		distNodes[1] = target;

	if(nof_distNodes > 2)
	{
		for (int i = 0; i < nof_nodes; ++i)
		{
			OrderedNodes[i].id = i;
			if(i == source || i == target)
				OrderedNodes[i].degree = -1;	
			else
				OrderedNodes[i].degree = Nodes[i+1] - Nodes[i];	
		}
		sort(OrderedNodes.begin(), OrderedNodes.end(), ONodesCompare);

		for (int i = 2; i < nof_distNodes; i++){
				distNodes[i] = OrderedNodes[i].id;
		}
	}

	return distNodes;
}




__device__ __forceinline__ int LaneID() {
    int ret;
    asm("mov.u32 %0, %laneid;" : "=r"(ret) );
    return ret;
}



__device__ int2 atomicMax(int2* address, int2 val)
{
    unsigned long long* addr_as_ull = (unsigned long long*)address;
    unsigned long long  old = *addr_as_ull;
    unsigned long long  assumed;
    do {
        assumed = old;
        int2* temp = (int2*)&assumed;
        if (val.x > temp->x || (val.x == temp->x && val.y > temp->y))
            old = atomicCAS(addr_as_ull, assumed, *(unsigned long long*)&val);
        else
            break;
    }while(assumed != old);
    return *((int2*)&old);
}


__device__ void atomicStore(int2* address, int2 val){
    unsigned long long* addr_as_ull = (unsigned long long*)address;
    unsigned long long  old = *addr_as_ull;
    unsigned long long  assumed;
    assumed = old;
    atomicCAS(addr_as_ull, assumed, *(unsigned long long*)&val);
    return;
}