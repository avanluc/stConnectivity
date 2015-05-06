#include "stConnectivity.hpp"

using namespace std;


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
			
			if(x >= N || y >= N)
				printf("Error at row %d: node id > # nodes\n", i+2);
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
	int *edges = (int*)calloc(E, sizeof(int));
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
				for(int l = (graph[i-1].x + 1); l <= N; l++)
					vertex[l] = vertex[graph[i-1].x] + j;
			j = 0;
		}
		if( i == ((E)-1) && graph[i].x < (N-1)){
			for (int l = (graph[i].x + 1); l <=N; l++)
				vertex[l]++;
		}
		// Fill adjacency list
		j++;
		edges[i] = graph[i].y;
		//printf("set arc %d->%d\n", graph[i].x, edges[i]);
	}

	// It is convenient to add a N-th element
	vertex[N] = E;
	printf("CSR structure created\n");
	printf("-------------------------\n");
	
	int source = atoi(argv[3]);
	int target = atoi(argv[4]);
	
	// Source/Target control
	if( source >= N )
	{
		printf("Source node > |N| = %d\n", N);
		return -1;
	}
	if( target >= N )
	{
		printf("Target node > |N| = %d\n", N);
		return -1;	
	}


	/*
	* St-Connectivity Algorithm 1
	*/
	if(argv[2][0] == '1')
	{
		for (int i = 0; i < 5; ++i)
		{
			clock_t start, end;
			int* Queue = (int*)calloc(N,sizeof(int));

			start = clock();

			bool connect = stConnectivity(vertex, edges, N, source, target, Queue);
			
			end = clock();
			// Calculate elapsed time
			double cpuTime = ((double)(end - start)) / CLOCKS_PER_SEC;

			printf("Elapsed time for st-Connectivity Algorithm 1 = %.3f s\n", cpuTime);
			printf("Result for st-Connectivity Algorithm 1 from %s to %s is %s\n", argv[3], argv[4], (connect ? "true" : "false"));
			free(Queue);
		}
	}
	/*
	* St-Connectivity Algorithm 2	
	*/
	else if(argv[2][0] == '2')
	{
		for (int i = 0; i < 5; ++i)
		{
			clock_t start, end;
			int* Queue = (int*)calloc(N,sizeof(int));
			int* Queue2 = (int*)calloc(N,sizeof(int));
			
			start = clock();
					
			//CALCOLO
			bool connect = BidirectionalStConnectivity(vertex, edges, N, source, target, Queue, Queue2);

			end = clock();
			// Calculate elapsed time
			double cpuTime = ((double)(end - start)) / CLOCKS_PER_SEC;
			printf("Elapsed time for t-Connectivity Algorithm 2 = %.3f s\n", cpuTime);
			printf("Result for st-Connectivity Algorithm 2 from %s to %s is %s\n", argv[3], argv[4], (connect ? "true" : "false"));
			free(Queue);
			free(Queue2);

		}
	}
	/*
	* St-Connectivity Algorithm 3
	*/
	else if(argv[2][0] == '3')
	{
		int nof_distNodes = atoi(argv[5]);
		for (int i = 0; i < 5; ++i)
		{
			bool* adjMatrix = (bool*)calloc(nof_distNodes*nof_distNodes, sizeof(bool));		//adjacence matrix distiguished nodes
			int** Queue = new int*[nof_distNodes];		// array of queues
			for(int i = 0; i < nof_distNodes; i++)
				Queue[i] = new int[N];

			clock_t start, end;
			start = clock();
					
			//CALCOLO
			bool connect = UlmannStConnectivity(vertex, edges, N, source, target, nof_distNodes, Queue, adjMatrix);

			end = clock();
			// Calculate elapsed time
			double cpuTime = ((double)(end - start)) / CLOCKS_PER_SEC;
			printf("Elapsed time for t-Connectivity Algorithm 3 = %.3f s\n", cpuTime);
			printf("Result for st-Connectivity Algorithm 3 from %s to %s is %s\n", argv[3], argv[4], (connect ? "true" : "false"));
			delete [] Queue;
			free(adjMatrix);
		}
	}
	else
	{
		int bfs = SimpleBFS(vertex, edges, N, source);
		printf("%d\n", bfs);
	}
	printf("-------------------------\n");
	// Free host memory
	free(vertex);
	free(edges);
	free(graph);
	return 0;
}