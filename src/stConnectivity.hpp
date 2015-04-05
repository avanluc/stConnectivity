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

/*
* Compare function for sorting
*/
int compare(const void *x, const void *y){
	edge a = *(edge *)x;
	edge b = *(edge *)y;
	return a.x < b.x ? -1 : (a.x > b.x ? 1 : (a.y < b.y ? -1 : a.y > b.y));
}


/*
* stConnectivity sequential function
*/
bool stConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
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
	delete Queue;
	return false;
}



/*
* Bidirectional StConnectivity sequential function 
*/
bool BidirectionalStConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int left2 = 0, right2 = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
	int* Queue2 = new int[nof_nodes];	// coda dei nodi da visitare
	//vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato
	vector<int> Visited(nof_nodes);	// indica se un nodo è stato visitato
	
	Queue[0] = source;					// Si parte dalla sorgente
	Queue2[0] = target;					// Si parte dalla sorgente
	Visited[source] = 1;				// si marca come visitata
	Visited[target] = 2;				// si marca come visitata
	
	while (left < right && left2 < right2) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
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
	delete Queue;
	return false;
}







/*
*  Ulmann's sequential version of stConnectivity function
*/
bool UlmannStConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target) {
	int nof_distNodes= (int) sqrt(nof_nodes)*log2(nof_nodes);
	int* randNodes = new int[nof_distNodes];
	// Initialize rand function
	srand (time(NULL));
	printf(" %d\n", nof_distNodes);

	// Controllare che i nodi siano distinti e aggiungere source e target 
	for(int i=1; i<nof_distNodes+1; i++)
	{
		randNodes[i-1] = (int) rand() % i;
		printf(" %d ", randNodes[i-1]);
	}


	// launch BFS from all distinguished nodes and return the list of distinguished nodes reachible from it.

	// create new graph H with only distinguished nodes and add arcs v->u if v can reach u.

	// Launch stConnectivity from source to target on H
	return false;
}



/*
* stConnectivity sequential function
*/
vector<int> BFS(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target, const int distance) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
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
	delete Queue;
	return false;
}