#pragma once

#include "stConn.h"


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
	std::vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato

	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata

	while (left < right)				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
	{
		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = 0; i < nof_nodes; ++i) 		// per ogni nodo nella lista di adiacenza di qNode
		{
			if(adjMatrix[qNode*nof_nodes + i])
				if (!Visited[i])					// se dest non è ancora stato visitato
				{ 
					Visited[i] = true;				// lo marco come visitato
					Queue[right++] = i;				// lo aggiungo alla coda dei nodi incontrati durante la visita 
				}
		}
	}
	return (Visited[target] == true);
}



/*
* Sort nodes per degree
*/
std::vector< ONodes > OrderNodes(const int* Nodes, const int nof_nodes){
	std::vector< ONodes > OrderedNodes(nof_nodes);
	for (int i = 0; i < nof_nodes; ++i)
	{
		OrderedNodes[i].id = i;
		OrderedNodes[i].degree = Nodes[i+1] - Nodes[i];	
	}
	sort(OrderedNodes.begin(), OrderedNodes.end(), ONodesCompare);
	return OrderedNodes;
}



/*
* Function that choose nof_distNodes nodes of the graph sorted by degree
*/
void ChooseNodes(int* sources, std::vector< ONodes > OrderedNodes, const int nof_distNodes, const int source, const int target) {

	sources[0] = source;
	
	if(nof_distNodes > 1)		
		sources[1] = target;
	if(nof_distNodes > 2)
		for (int i = 2; i < nof_distNodes; i++)
			sources[i] = OrderedNodes[i].id;
	return;
}



/*
* Function that choose nof_distNodes nodes of the graph
*/
void ChooseRandomNodes(int* sources, const int* Nodes, const int nof_nodes, const int nof_distNodes, const int source, const int target) {

	sources[0] = source;
	
	if(nof_distNodes > 1)
		sources[1] = target;

	if(nof_distNodes > 2)
	{
		for (int i = 2; i < nof_distNodes; i++)
		{
			bool isGood = 1;
			sources[i] = rand() % nof_nodes;
			for(int j = 0; j < i; j++)
			{
				if(sources[j] == sources[i])
				{
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



/*
* Function to convert graph into CSR structure
* !!! DEPRECATED !!!
*/
void GraphToCSR(const edge* graph, int* vertex, int* edges, int N, int E){
	int j = 0;

	for(int i = 0; i < (E); i++)
	{
		if(i!=0 && graph[i].x != graph[i-1].x)
		{
			if( (graph[i].x - graph[i-1].x) == 1 )
				vertex[graph[i].x] = vertex[graph[i-1].x] + j;
			else
				for(int l = (graph[i-1].x + 1); l <= N; l++)
					vertex[l] = vertex[graph[i-1].x] + j;
			j = 0;
		}
		if( i == ((E)-1) && graph[i].x < (N-1))
			for (int l = (graph[i].x + 1); l <=N; l++)
				vertex[l]++;
		j++;
		edges[i] = graph[i].y;
	}
	vertex[N] = E;
	return;
}



/*
* Function to read graph from file
* !!! DEPRECATED !!!
*/
void ReadGraph(char *file, edge *graph, int *N, int *E){
	int x,y;
	std::ifstream in (file);
	if(in.is_open()){
		in >> *N >> *E;
		for(int i = 0; i < *E; i++)
		{	
			in >> x >> y;
			graph[i].x = x;
			graph[i].y = y;
			if(x >= *N || y >= *N)
				printf("Error at row %d: node id > # nodes\n", i+2);
		}
		// Sorting graph using specified compare function
		qsort(graph, (*E), sizeof(edge), compare);
	}
	return;
}