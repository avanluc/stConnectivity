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
	delete Queue;
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









bool UlmannStConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target, const int nof_distNodes) {

	bool* adjMatrix = (bool*)calloc(nof_distNodes*nof_distNodes, sizeof(bool));		//adjacence matrix distiguished nodes
	vector< int > distNodes(nof_distNodes);			// vector of distinguished nodes
	vector< ONodes > OrderedNodes(nof_nodes);		// vector of nodes ordered by degree
	//srand (time(NULL));			// Initialize rand function

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

	int** Queue = new int*[nof_distNodes];		// array of queues
	vector<int> Marker(nof_nodes);				// vector of markers
	vector<int> left(nof_distNodes);
	vector<int> right(nof_distNodes);
	int finishedBFS = 0;

	// initialization loop
	for(int i = 0; i < nof_distNodes; i++)
	{
	    Queue[i] = new int[nof_nodes];
	    Queue[i][0] = distNodes[i];
	    left[i] = 0; right[i] = 1;
	    Marker[distNodes[i]] = distNodes[i]+1;	// every distNode marks with (distNodeId + 1)
	}

	while( finishedBFS < nof_distNodes )
	{
		//printf("\n--New Step For All--\n\n");
		// nof_distNodes bfs that advance step by step
		for (int i = 0; i < nof_distNodes; i++)
		{
			if(left[i] < right[i])
			{
				int qNode = Queue[i][left[i]++];		// qNode = prossimo nodo da cui far partire la visita
				//printf("Executing a BFS step from %d\n", distNodes[i]);
				//printf("looking for %d neigbors\n", qNode);
				for (int j = Nodes[qNode]; j < Nodes[qNode + 1]; j++) // per ogni nodo nella lista di adiacenza di qNode
				{
					int dest = Edges[j];										// dest = destinazione arco
					if (Marker[dest] == 0) {						// se dest non è ancora stato visitato 
						Marker[dest] = distNodes[i]+1;				// lo marco come visitato
						Queue[i][right[i]++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita  
						//printf("\tvisited %d\n", dest);
					}
					else if(Marker[dest] != distNodes[i]+1){
						// aggiungere Marker[dest] alla matrice!!!!
						int markerPos = find(distNodes.begin(), distNodes.end(), Marker[dest]-1) - distNodes.begin();
						adjMatrix[markerPos + nof_distNodes*i] = true;
						adjMatrix[i + nof_distNodes*markerPos] = true;
						
						//printf("\tnode %d already visited by %d\n", dest, Marker[dest]-1);
						// printf("\tadding arc between %d and %d\n", i, markerPos);
						// printf("\tsetting true to adjMatrix[%d] and adjMatrix[%d]\n", markerPos + nof_distNodes*i, i + nof_distNodes*markerPos);
					}
				}	
			}
			else if (left[i] == right[i]){
				finishedBFS++; left[i]++;
				//printf("Node %d has finished!!\n", distNodes[i]);
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

	free(adjMatrix);
	delete [] Queue;
	return b;
}








/*
* stConnectivity sequential function
*/
/*
vector<int> BFS(const int* Nodes, const int* Edges, const int nof_nodes,  int source,  int* isRand,  int distance) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	int* Queue = new int[nof_nodes];	// coda dei nodi da visitare
	vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato
	vector<int> distances(nof_nodes);
	vector<int> adj;
	
	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata
	distances[source] = 0;

	//printf("Executing BFS from %d\n", source);
	int d = 0;

	while (left < right && d < distance) {				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
		//printf("IN while\n");

		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = Nodes[qNode]; i < Nodes[qNode + 1]; ++i) {		// per ogni nodo nella lista di adiacenza di qNode
			int dest = Edges[i];										// dest = destinazione arco

			if (!Visited[dest]) {					// se dest non è ancora stato visitato 
				Visited[dest] = true;				// lo marco come visitato
				distances[dest] = distances[qNode] +1;
				Queue[right++] = dest;				// lo aggiungo alla coda dei nodi incontrati durante la visita 
				if ( distances[dest] > distance )
					d = distance +1;
				//printf("\tvisited %d at distance %d\n", dest, distances[dest]);
				//printf("\tarc %d --> %d at distance %d\n",qNode, dest, distances[dest]);
			}
		}
	}
	for (int i = 0; i < nof_nodes; i++)
		if( Visited[i] && i != source && isRand[i])
			adj.push_back(i);
	
	delete Queue;
	return adj;
	//return false;
}

*/





/*
*  Ulmann's sequential version of stConnectivity function
*/
/*
bool UlmannStConnectivity1(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target) {

	//int nof_distNodes= (int) sqrt(nof_nodes)*log2(nof_nodes);
	int nof_distNodes= (int) sqrt(nof_nodes);
	int* isRand = new int[nof_nodes];
	int count = nof_distNodes;
	int distEdges = 0;

	vector< int > distNodes(nof_distNodes);
	//vector< vector<int> > adjLists(nof_nodes);
	vector< vector<int> > adjLists(nof_distNodes);
	//vector<int> utility;
	srand (time(NULL));			// Initialize rand function


	distNodes[0] = source;
	distNodes[1] = target;
	isRand[source] = 1;
	isRand[target] = 1;

	for (int i = 2; i < nof_distNodes; i++)
	{
		int n = (int)rand() % nof_nodes;
		if(find(distNodes.begin(), distNodes.end(), n)==distNodes.end()){
			distNodes[i] = n;
			isRand[n] = 1;
		}
	}

	// Controllare che i nodi siano distinti e aggiungere source e target 
	// for (int i = 0; i < nof_nodes; i++)
	// {
	// 	isRand[i] = (int) rand() % 2;
	// 	if(isRand[i] == 1)
	// 		count++;
	// }
	// if (!isRand[source])
	// {
	// 	isRand[source] = 1;
	// 	count++;
	// }
	// if (!isRand[target])
	// {
	// 	nof_distNodes
	// 	count++;
	// }

	// printf("N = %d,  count = %d\n",nof_nodes, count);
	
	// launch BFS from all distinguished nodes and return the list of distinguished nodes reachible from it.
	// for (int i = 0; i < nof_nodes; i++)
	// {
	// 	if (isRand[i])
	// 	{
	// 		adjLists[i] = BFS(Nodes, Edges, nof_nodes, i, isRand, count);
	// 		distEdges += adjLists[i].size();
	// 		//utility.push_back(i);
	// 	}
	// 	else
	// 		adjLists[i].push_back(-1);
	// }
	// for (int i = 0; i < nof_distNodes; i++)
	// {
	// 	adjLists[i] = BFS(Nodes, Edges, nof_nodes, distNodes[i], isRand, count);
	// 	distEdges += adjLists[i].size();
	// }
	// printf("N = %d,  count = %d\n",nof_nodes, count);
	
	// create new graph H with only distinguished nodes and add arcs v->u if v can reach u.
	int *H_vertex = (int*)calloc(count, sizeof(int));
	int *H_edges = (int*)calloc(distEdges, sizeof(int));

	int k = 0;
	int new_source = 0; 
	int new_target = 1;

	for(int i = 0; i < nof_distNodes; i++)
	{
		H_vertex[i] = k;
		for (int j = 0; j < adjLists[i].size(); j++)
			H_edges[k++] = find(distNodes.begin(), distNodes.end(), adjLists[i][j]) - distNodes.begin();
	}

	// int k = 0; int l = 0;
	// int new_source = 0; 
	// int new_target = 1;

	// for(int i = 0; i < nof_nodes; i++)
	// {
	// 	if(isRand[i])
	// 	{
	// 		// if (i == source)
	// 		// 	new_source = l;
	// 		// if (i == target)
	// 		// 	new_target = l;
	// 		H_vertex[l++] = k;
	// 		for (int j = 0; j < adjLists[i].size(); j++)
	// 		{
	// 			//H_edges[k++] = adjLists[i][j];
	// 			H_edges[k++] = find(distNodes.begin(), distNodes.end(), adjLists[i][j]) - distNodes.begin();
	// 		}
	// 	}
	// }
	// printf("N = %d,  count = %d\n",nof_nodes, count);
	// printf("\n\n");
	// for (int i = 0; i < count; i++)
	// 	printf("| %d ", H_vertex[i]);
	// printf("\n\n");
	// for (int i = 0; i < distEdges; i++)
	// 	printf("| %d ", H_edges[i]);
	// printf("\n\n");
	
	// Launch stConnectivity from source to target on H
	bool b = stConnectivity(H_vertex, H_edges, count, new_source, new_target);
	free(H_vertex);
	free(H_edges);
	free(isRand);
	return b;
}

*/

