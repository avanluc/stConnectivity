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
* Min function for percentage evaluation
*/
double min(std::vector<long double> data, int n){
	double min  = 100;
	for(int i = 0; i < n; i++){
		min  = (data[i] < min ? data[i] : min);
	}
	return min;
}



/*
* Max function for percentage evaluation
*/
double max(std::vector<long double> data, int n){
	double max  = 0;
	for(int i = 0; i < n; i++){
		max  = (data[i] > max ? data[i] : max);
	}
	return max;
}



/*
* Evaluate mean percentage of visited graph
*/
void computeMeanPercentage(std::vector<long double> Percentual, int percentCnt){
	double sum = 0;
	for(int i = 0; i < percentCnt; i++)
		sum += Percentual[i];

	printf("# Completed Visit: %d on %d\n", (N_TEST - percentCnt), N_TEST);
	if(percentCnt != 0)
	{
		printf("AVG Percentual \t\t: %.2f%\n", sum / percentCnt);
		printf("MIN Percentual \t\t: %.2f%\n", min(Percentual, percentCnt));
		printf("MAX Percentual \t\t: %.2f%\n", max(Percentual, percentCnt));
	}
}



/*
* Evaluate Elapsed Time
*/
void computeElapsedTime(std::vector<double> par_times, std::vector<double> seq_times, int connectCnt){
	double sum_par = 0;
	double sum_seq = 0;

	for (int i = 1; i < N_TEST; ++i){
		sum_par += par_times[i];
		sum_seq += seq_times[i];
	}

	printf("\n# Positive Responses: %d on %d\n", (N_TEST - connectCnt), N_TEST);
	printf("AVG TIME \t\t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, (sum_par + sum_seq) / (N_TEST), 27, 0);
	printf("AVG PARALLEL TIME \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_par / (N_TEST), 27, 0);
	printf("AVG MATRIX BFS TIME \t: %c[%d;%dm%.1f%c[%dm ms\n\n", 27, 0, 31, sum_seq / (N_TEST), 27, 0);
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