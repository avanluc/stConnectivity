#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <utility>
#include <climits>
#include <vector>

using namespace std;

// Graph parameters (#nodes, #edges)
int  N, E;


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



// Grafo passato in forma CRS
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



// Grafo passato in forma CRS
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





int main(int argc, char *argv[]){

	// If wrong number of argumants print usage
	if(argc < 3)
	{
		printf("Usage: ./stConnectivity 'input_file' '1/2/3' 'source' 'target'\n");
		return -1;
	}


	// Read graph from file
	int x,y;
	edge *graph; 
	ifstream in (argv[1]);
	if(in.is_open())
	{
		// Read graph parameters and allocate memory
		in >> N >> E;
		graph = (edge*)malloc(E*sizeof(edge));

		// Read graph
		for(int i = 0; i < E; i++)
		{	
			in >> x >> y;
			graph[i].x = x;
			graph[i].y = y;

			// // Duplicate the arc
			// graph[i+E].x = y;
			// graph[i+E].y = x;
		}

		// Sorting graph using specified compare function
		qsort(graph, (E), sizeof(edge), compare);
	}
	else
	{
		printf("No input file\n");
		return -1;
	}
	in.close();
	printf("Graph sorted\n");

	// for(int i=0; i<20; i++)
	// 	printf("%d, %d\n", graph[i].x, graph[i].y);

	// Memory allocation and initialization
	int *vertex = (int*)calloc(N+1, sizeof(int));
	int *edges = (int*)calloc(E*2, sizeof(int));
	printf("Memory allocated\n");

	
	// Creation of CSR structure

	int j = 0;  // Current adjacency list length

	// for each arc i in the graph
	for(int i = 0; i < (E); i++)
	{
		// if i isn't the first arc and 
		// the src of i is different from the src of i-1 
		if(i!=0 && graph[i].x != graph[i-1].x)
		{
			// if the difference between arcs sources == 1 then set src pointer
			if( (graph[i].x - graph[i-1].x) == 1 )
				vertex[graph[i].x] = vertex[graph[i-1].x] + j;
			else
				// else set each pointer between the sources to j
				for(int l = (graph[i-1].x + 1); l <= graph[i].x; l++)
					vertex[l] = vertex[graph[i-1].x] + j;
			j = 0;
		}
		// Fill adjacency list
		j++;
		// edges[i].x = graph[i].x;
		// edges[i].y = graph[i].y;
		edges[i] = graph[i].y;
	}

	// It is convenient to add a N-th element
	vertex[N] = (E);
	printf("CSR structure created\n");
	
	// Source/Target control
	if( atoi(argv[3]) >= N )
	{
		printf("Source node > |N| = %d\n", N);
		return -1;
	}
	if( atoi(argv[4]) >= N )
	{
		printf("Target node > |N| = %d\n", N);
		return -1;	
	}

	/*
	* St-Connectivity Algorithm 1
	*/
	if(argv[2][0] == '1')
	{
		clock_t start, end;
		start = clock();
		bool connect = stConnectivity(vertex, edges, N, atoi(argv[3]), atoi(argv[4]));
		
		end = clock();
		// Calculate elapsed time
		double cpuTime = ((double)(end - start)) / CLOCKS_PER_SEC;

		printf("\nElapsed time for st-Connectivity Algorithm 1 = %.3f s\n", cpuTime);
		printf("Result for st-Connectivity Algorithm 1 from %s to %s is %s\n", argv[3], argv[4], (connect ? "true" : "false"));

	}
	/*
	* St-Connectivity Algorithm 2
	*/
	else if(argv[2][0] == '2')
	{
		clock_t start, end;
		start = clock();
				
		//CALCOLO
		bool connect = BidirectionalStConnectivity(vertex, edges, N, atoi(argv[3]), atoi(argv[4]));

		end = clock();
		// Calculate elapsed time
		double cpuTime = ((double)(end - start)) / CLOCKS_PER_SEC;

		printf("Elapsed time for t-Connectivity Algorithm 2 = %.3f s\n", cpuTime);
		printf("Result for st-Connectivity Algorithm 2 from %s to %s is %s\n", argv[3], argv[4], (connect ? "true" : "false"));
		

		// Print results
	}
	/*
	* St-Connectivity Algorithm 3
	*/
	else
	{
		clock_t start, end;
		start = clock();
				
		//CALCOLO


		end = clock();
		// Calculate elapsed time
		double cpuTime = ((double)(end - start)) / CLOCKS_PER_SEC;

		printf("Elapsed time for t-Connectivity Algorithm 3 = %.3f s\n", cpuTime);
		

		// Print results
	}

	// Free host memory
	free(vertex);
	free(edges);
	free(graph);
	
	return 0;
}