#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include "definition.cuh"

using namespace std;

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


/*
* BFS on adjacency matrix performed on CPU
*/
bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target, int* Queue) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	//int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
	vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato
	
	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata

	while (left < right) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = 0; i < nof_nodes; ++i) 		// per ogni nodo nella lista di adiacenza di qNode
			if(adjMatrix[qNode*nof_nodes + i])
				if (!Visited[i]) {					// se dest non è ancora stato visitato 
					Visited[i] = true;				// lo marco come visitato
					Queue[right++] = i;				// lo aggiungo alla coda dei nodi incontrati durante la visita 
				}
	}
	//delete [] Queue;
	return (Visited[target] == true);
}


/*
* Function that choose nof_distNodes nodes of the graph
*/
void ChooseNodes(int* sources, const int* Nodes, const int nof_nodes, const int nof_distNodes, const int source, const int target) {

	//vector< int > distNodes(nof_distNodes);			// vector of distinguished nodes
	//vector< ONodes > OrderedNodes(nof_nodes);		// vector of nodes ordered by degree

	sources[0] = source;
	
	if(nof_distNodes > 1)
		sources[1] = target;

	if(nof_distNodes > 2)
	{
		/*for (int i = 0; i < nof_nodes; ++i)
		{
			OrderedNodes[i].id = i;
			if(i == source || i == target)
				OrderedNodes[i].degree = -1;	
			else
				OrderedNodes[i].degree = Nodes[i+1] - Nodes[i];	
		}
		sort(OrderedNodes.begin(), OrderedNodes.end(), ONodesCompare);*/
		//srand (time(NULL));
		for (int i = 2; i < nof_distNodes; i++){
			bool isGood = 1;
			sources[i] = rand() % N;
			for(int j = 0; j < i; j++){
				if(sources[j] == sources[i]){
					isGood = 0;
					break;
				}
			}
			if (!isGood)
				i--;
		}
	}
	return;
}


void GraphToCSR(const edge* graph, int* vertex, int* edges){
	int j = 0;  // Current adjacency list length

	// for each arc i in the graph
	for(int i = 0; i < (E); i++){
		// if i isn't the first arc and 
		// the src of i is different from the src of i-1 
		if(i!=0 && graph[i].x != graph[i-1].x){
			// if the difference between arcs sources == 1 then set src pointer
			if( (graph[i].x - graph[i-1].x) == 1 )
				vertex[graph[i].x] = vertex[graph[i-1].x] + j;
			// else set each pointer between the sources to j
			else
				for(int l = (graph[i-1].x + 1); l <= N; l++)
					vertex[l] = vertex[graph[i-1].x] + j;
			j = 0;
		}
		if( i == ((E)-1) && graph[i].x < (N-1))
			for (int l = (graph[i].x + 1); l <=N; l++)
				vertex[l]++;
		// Fill adjacency list
		j++;
		edges[i] = graph[i].y;
	}
	// It is convenient to add a N-th element
	vertex[N] = E;
	return;
}


void ReadGraph(char *file, edge *graph){
	int x,y;
	ifstream in (file);
	if(in.is_open()){
		in >> N >> E;
		for(int i = 0; i < E; i++){	
			in >> x >> y;
			graph[i].x = x;
			graph[i].y = y;
			if(x >= N || y >= N)
				printf("Error at row %d: node id > # nodes\n", i+2);
		}
		// Sorting graph using specified compare function
		qsort(graph, (E), sizeof(edge), compare);
	}
	return;
}