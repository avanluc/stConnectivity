#include <KernelStConnectivity.hpp>

using namespace std;


int main(int argc, char *argv[]){

	// If wrong number of argumants print usage
	if(argc < 3)
	{
		printf("Usage: ./stConnectivity 'input_file' 'source' 'target' '#_distinguished_nodes'\n");
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
	//  	printf("%d, %d\n", graph[i].x, graph[i].y);

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

	int source = atoi(argv[2]);
	int target = atoi(argv[3]);
	int nof_distNodes = atoi(argv[4]);
	
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
	
	// calculate size for allocation
	size_t sizeN = (N+1) * sizeof(int);
	size_t sizeE = E * sizeof(int);
	size_t sizeN2 = N * sizeof(int2);
	size_t sizeN_dist = nof_distNodes * nof_distNodes * sizeof(bool);

	// Allocate adj matrix in row major
	bool *matrix = (bool*)calloc(nof_distNodes * nof_distNodes, sizeof(bool));
	int2 *Dist_Col = (int2*)calloc(N, sizeof(int2));
	bool *newlevel = (bool*)calloc(1, sizeof(bool));
	
	// Calculate Kernel grid dimension
	int dimGrid = ((N)/BLOCK_SIZE)+1;

	// Allocate device memory
	int *Dedges;
	int *Dvertex;
	bool *DMatrix;
	int2 *Ddist_Col;
	bool *DnewLevel;

	gpuErrchk( cudaMalloc((void **) &Dvertex, sizeN) );
	gpuErrchk( cudaMalloc((void **) &Dedges, sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, sizeN_dist) );
	gpuErrchk( cudaMalloc((void **) &Ddist_Col, sizeN2) );
	gpuErrchk( cudaMalloc((void **) &DnewLevel, sizeof(bool)) );
	printf("Device memory allocated\n");
 
    // choose nof_distNodes distinguished nodes with source and target in it and return a vector of them	
    vector< int > id = ChooseNodes(vertex, N, nof_distNodes, source, target);		// !!!!!!!!!!
    printf("Distinguished nodes chosen\n");
    // inizializzazione dei valori a INT_MAX
    for (int i = 0; i < N; i++)
    {
    	Dist_Col[i].x = INT_MAX;
		Dist_Col[i].y = INT_MAX;
    }

    // inizializzazione dei valori dei nodi distinti a Distance 0 e Color id
    for (int i = 0; i < nof_distNodes; i++)
	{
		//printf("distinguished node %d = %d\n", i, id[i]);
		int j = id[i];
		Dist_Col[j].x = 0;
		Dist_Col[j].y = i;
	}
	printf("Valori inizializzati\n");

    // Copy host memory for vertex, edges and results vectors to device
    gpuErrchk( cudaMemcpy(Dvertex, vertex, sizeN, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(Dedges, edges, sizeE, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(DMatrix, matrix, sizeN_dist, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(Ddist_Col, Dist_Col, sizeN2, cudaMemcpyHostToDevice) );
    printf("Memcopy executed\n");

	// current level of visit
    int level = 0;
    newlevel[0] = true;

	// Allocate CUDA events to be used for timing
    cudaEvent_t start;
    cudaEvent_t stop;

    gpuErrchk( cudaEventCreate(&start) );
    gpuErrchk( cudaEventCreate(&stop) );
    
    // Record the start event
    gpuErrchk( cudaEventRecord(start, NULL) );

    // Launch Cuda Kernel
	dim3 block(BLOCK_SIZE, 1);
    dim3 grid(dimGrid, 1);

    // Problema nel ritornare il valore di newlevel!!!!! 
    // dovrebbe essere una variabile sul device passata come parametro e quindi ritornate tramite memcopy

    //printf("Begin visit with cuda Kernels\n");
    while(newlevel[0])
    {
    	newlevel[0] = false;
    	
    	gpuErrchk( cudaMemcpy(DnewLevel, newlevel, sizeof(bool), cudaMemcpyHostToDevice) );
    	
    	//printf("\nCalling kernel with parameters(Dvertex, Dedges, Ddist_Col, %d, false, DMatrix, %d, %d)", level, N, nof_distNodes);
    	stConn<<< grid, block >>>(Dvertex, Dedges, Ddist_Col, level, DnewLevel, DMatrix, N, nof_distNodes);
    	//printf("\nExecuted a visit to level %d\n\n", level);
		
		gpuErrchk( cudaMemcpy(newlevel, DnewLevel, sizeof(bool), cudaMemcpyDeviceToHost) );
    	
    	//printf("newlevel = %d\n", newlevel[0]);
    	level++;
    }

    // Provare a vedere se è possibile lanciare un secondo kernel per calcolare stConn sulla matrice 
    // mantenendo quest'ultima in memoria device e dopo ritornare solo un booleano
    
    // Copy result vector from device to host
	gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeN_dist, cudaMemcpyDeviceToHost) );


	// printf("matrix completed\n");

	// // Print matrix
	// for (int i = 0; i < nof_distNodes; ++i)
	// {
	// 	printf("| ");
	// 	for (int j = 0; j < nof_distNodes; ++j)
	// 		printf("%d ", matrix[nof_distNodes*i+j]);
	// 		//printf("i:%d, j:%d, position: %d \n", i, j, nof_distNodes*i+j);
	// 	printf("|\n");
	// }
	//printf("Calling MatrixBFS with arguments(matrix, %d, %d, %d)\n", nof_distNodes, 0, 1);

	// altrimenti, calculate stConn on the adj matrix on host
	bool connect = MatrixBFS(matrix, nof_distNodes, 0, 1);

    //return;
	gpuErrchk( cudaEventRecord(stop, NULL) );
    // Wait for the stop event to complete
    gpuErrchk( cudaEventSynchronize(stop) );

    // Claculate elapsed time
    float msecTotal = 0.0f;
    gpuErrchk( cudaEventElapsedTime(&msecTotal, start, stop) );
    msecTotal /= 1000;

	printf("Elapsed time for t-Connectivity Algorithm 3 = %.3f s\n", msecTotal);
	printf("Result for st-Connectivity Algorithm 3 from %d to %d is %s\n", source, target, (connect ? "true" : "false"));
	printf("-------------------------\n");
	
	//free memoria device
	cudaFree(Dvertex);
    cudaFree(Dedges);
    cudaFree(DMatrix);
    cudaFree(Ddist_Col);
    cudaFree(DnewLevel);

	free(matrix);
	free(Dist_Col);
	free(vertex);
	free(edges);
	free(graph);
	free(newlevel);

	return 0;
}