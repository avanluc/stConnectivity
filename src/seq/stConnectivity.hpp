#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <utility>
#include <climits>
#include <vector>
#include <math.h> 

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


int SimpleBFS(const int* Nodes, const int* Edges, const int nof_nodes, const int source) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
	vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato
	
	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata
	
	while (left < right) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = Nodes[qNode]; i < Nodes[qNode + 1]; ++i) {		// per ogni nodo nella lista di adiacenza di qNode
			int dest = Edges[i];										// dest = destinazione arco
			if (!Visited[dest]) {					// se dest non è ancora stato visitato 
				Visited[dest] = true;				// lo marco come visitato
				Queue[right++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita  
			}
		}
	}
	return Queue[right-1];
}


/*
* stConnectivity sequential function
*/
bool stConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target, int* Queue) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato
	
	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata
	
	while (left < right) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = Nodes[qNode]; i < Nodes[qNode + 1]; ++i) {		// per ogni nodo nella lista di adiacenza di qNode
			int dest = Edges[i];										// dest = destinazione arco

			if (dest == target)						// se dest == target return true
				return true;
			if (!Visited[dest]) {					// se dest non è ancora stato visitato 
				Visited[dest] = true;				// lo marco come visitato
				Queue[right++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita  
			}
		}
	}
	return false;
}



/*
* Bidirectional StConnectivity sequential function 
*/
bool BidirectionalStConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target, int* Queue, int* Queue2) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int left2 = 0, right2 = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	vector<int> Visited(nof_nodes);	// indica se un nodo è stato visitato
	
	Queue[0] = source;					// Si parte dalla sorgente
	Queue2[0] = target;					// Si parte dalla sorgente
	Visited[source] = 1;				// si marca come visitata
	Visited[target] = 2;				// si marca come visitata
	
	while (left < right && left2 < right2) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
		if(right <= right2)
		{
			int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

			for (int i = Nodes[qNode]; i < Nodes[qNode + 1]; ++i) {		// per ogni nodo nella lista di adiacenza di qNode
				int dest = Edges[i];										// dest = destinazione arco

				if (dest == target || Visited[dest] == 2)						// se dest == target return true
					return true;
				if (Visited[dest] == 0) {					// se dest non è ancora stato visitato 
					Visited[dest] = 1;				// lo marco come visitato
					Queue[right++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita  
				}
			}
		}
		else
		{
			int qNode2 = Queue2[left2++];		// qNode = prossimo nodo da cui far partire la visita

			for (int i = Nodes[qNode2]; i < Nodes[qNode2 + 1]; ++i) {		// per ogni nodo nella lista di adiacenza di qNode
				int dest = Edges[i];										// dest = destinazione arco

				if (dest == source || Visited[dest] == 1)						// se dest == target return true
					return true;
				if (Visited[dest] == 0) {					// se dest non è ancora stato visitato 
					Visited[dest] = 2;				// lo marco come visitato
					Queue2[right2++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita  
				}
			}
		}
	}
	return false;
}





/*
* stConnectivity sequential function
*/

bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
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
	
	delete Queue;
	return Visited[target];
}









bool UlmannStConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target, const int nof_distNodes, int ** Queue, bool* adjMatrix) {

	
	vector< int > distNodes(nof_distNodes);			// vector of distinguished nodes
	vector< ONodes > OrderedNodes(nof_nodes);		// vector of nodes ordered by degree

	distNodes[0] = source;
	distNodes[1] = target;

	for (int i = 0; i < nof_nodes; ++i)
	{
		OrderedNodes[i].id = i;
		OrderedNodes[i].degree = Nodes[i+1] - Nodes[i];
	}
	sort(OrderedNodes.begin(), OrderedNodes.end(), ONodesCompare);

	for (int i = 2; i < nof_distNodes; i++){
		if(OrderedNodes[i].id != source && OrderedNodes[i].id != target)
			distNodes[i] = OrderedNodes[i].id;
		else
			distNodes[i] = OrderedNodes[i+nof_distNodes].id;
	}

	//printf("Random nodes has been chosen\n");

	vector<int> Marker(nof_nodes);				// vector of markers
	vector<int> left(nof_distNodes);
	vector<int> right(nof_distNodes);
	int finishedBFS = 0;

	// initialization loop
	for(int i = 0; i < nof_distNodes; i++)
	{
	    Queue[i][0] = distNodes[i];
	    left[i] = 0; right[i] = 1;
	    Marker[distNodes[i]] = distNodes[i]+1;	// every distNode marks with (distNodeId + 1)
	}

	while( finishedBFS < nof_distNodes )
	{
		// nof_distNodes bfs that advance step by step
		for (int i = 0; i < nof_distNodes; i++)
		{
			if(left[i] < right[i])
			{
				int qNode = Queue[i][left[i]++];		// qNode = prossimo nodo da cui far partire la visita
				for (int j = Nodes[qNode]; j < Nodes[qNode + 1]; j++) // per ogni nodo nella lista di adiacenza di qNode
				{
					int dest = Edges[j];										// dest = destinazione arco
					if (Marker[dest] == 0) {						// se dest non è ancora stato visitato 
						Marker[dest] = distNodes[i]+1;				// lo marco come visitato
						Queue[i][right[i]++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita  
					}
					else if(Marker[dest] != distNodes[i]+1){
						int markerPos = find(distNodes.begin(), distNodes.end(), Marker[dest]-1) - distNodes.begin();
						adjMatrix[markerPos + nof_distNodes*i] = true;
						adjMatrix[i + nof_distNodes*markerPos] = true;
					}
				}	
			}
			else if (left[i] == right[i]){
				finishedBFS++; left[i]++;
			}
		}
	}	

	//printf("matrix completed\n");

	// // Print matrix
	// for (int i = 0; i < nof_distNodes; ++i)
	// {
	// 	printf("| ");
	// 	for (int j = 0; j < nof_distNodes; ++j)
	// 		printf("%d ", adjMatrix[nof_distNodes*i+j]);
	// 		//printf("i:%d, j:%d, position: %d \n", i, j, nof_distNodes*i+j);
	// 	printf("|\n");
	// }
	
	// // Launch stConnectivity from source to target on adjmatrix
	bool b = MatrixBFS(adjMatrix, nof_distNodes, 0, 1);

	return b;
}
