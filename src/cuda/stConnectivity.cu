#include <KernelStConnectivity_WEfficient.cuh>

using namespace std;


int main(int argc, char *argv[]){

	// If wrong number of argumants print usage
	if(argc < 3)
	{
		printf("\nUsage: ./stConnectivity 'input_file' '#_distinguished_nodes'\n\n");
		return -1;
	}

	// Read graph from file
	edge *graph;
	ifstream in (argv[1]);
	in >> N >> E;
	graph = (edge*)malloc(E*sizeof(edge));
	ReadGraph(argv[1], graph);
	
	// Memory allocation and initialization
	int *vertex = (int*)calloc(N+1, sizeof(int));
	int *edges = (int*)calloc(E, sizeof(int));
	
	// Creation of CSR structure
	GraphToCSR(graph, vertex, edges);
	//printf("Transfert graph to CSR structure\n");
	
	// int source = atoi(argv[2]);
	// int target = atoi(argv[3]);
	// int nof_distNodes = atoi(argv[4]);
	
	// // Source/Target control
	// if( source >= N )
	// {
	// 	printf("Source node > |N| = %d\n", N);
	// 	return -1;
	// }
	// if( target >= N )
	// {
	// 	printf("Target node > |N| = %d\n", N);
	// 	return -1;	
	// }

	int nof_distNodes = atoi(argv[2]);

	// calculate size for allocation
	size_t sizeN 	= (N+1) * sizeof(int);
	size_t sizeE 	= E * sizeof(int);
	size_t sizeN2 	= N * sizeof(int2);
	size_t sizeN3 	= nof_distNodes * sizeof(int);
	size_t sizeMATRIX = nof_distNodes * nof_distNodes * sizeof(bool);

	// Allocate adj matrix in row major
	bool *matrix 	= (bool*)calloc(nof_distNodes * nof_distNodes, sizeof(bool));
	int2 *Dist_Col 	= (int2*)calloc(N, sizeof(int2));
	int *Distance 	= (int*)calloc(nof_distNodes, sizeof(int));	
	int *sources 	= (int*)calloc(nof_distNodes, sizeof(int));
	int *Queue 		= (int*)calloc(nof_distNodes, sizeof(int));	

	// Allocate device memory
	int *Dedges;
	int *Dvertex;
	int *Ddistance;
	int *Dsources;
	int2 *Ddist_Col;
	bool *DMatrix;
	//bool *DnewLevel;

	gpuErrchk( cudaMalloc((void **) &Dvertex, sizeN) );
	gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN3) );
	gpuErrchk( cudaMalloc((void **) &Dedges, sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, sizeMATRIX) );
	gpuErrchk( cudaMalloc((void **) &Ddist_Col, sizeN2) );
	gpuErrchk( cudaMalloc((void **) &Dsources, sizeN3) );
	//printf("Device memory allocated\n");
	//printf("-----------------------------------\n");

 	vector<double> mean_times(3);
 	vector<double> par_times(N_TEST);
 	vector<double> seq_times(N_TEST);

	srand (time(NULL));
 	for (int test = 0; test < N_TEST; test++)
	{

		int source = rand() % N;
		int target = rand() % N;

		// choose nof_distNodes distinguished nodes with source and target in it and return a vector of them	
	    ChooseNodes(sources, vertex, N, nof_distNodes, source, target);

	    // inizializzazione dei valori a INT_MAX
	    for (int i = 0; i < N; i++)
	    {
	    	Dist_Col[i].x = INT_MAX;
			Dist_Col[i].y = INT_MAX;
	    }

	    // inizializzazione dei valori dei nodi distinti a Distance 0 e Color id
	    for (int i = 0; i < nof_distNodes; i++)
		{
			int j = sources[i];
			Dist_Col[j].x = 0;
			Dist_Col[j].y = i;
			Distance[i] = INT_MAX;
		}
		Distance[0] = 0;
		//bool NextLevel[4] = {0,0,0,0};
		
	    // Copy host memory for vertex, edges and results vectors to device
	    gpuErrchk( cudaMemcpy(Dvertex, vertex, sizeN, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Dedges, edges, sizeE, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(DMatrix, matrix, sizeMATRIX, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Ddist_Col, Dist_Col, sizeN2, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Ddistance, Distance, sizeN3, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Dsources, sources, sizeN3, cudaMemcpyHostToDevice) );
	    //gpuErrchk( cudaMemcpyToSymbolAsync(devNextLevel, NextLevel, sizeof(bool)*4, 0,cudaMemcpyHostToDevice) );

		// Allocate CUDA events to be used for timing
	    cudaEvent_t start;
	    cudaEvent_t start1;
	    cudaEvent_t stop;
	    cudaEvent_t stop1;

	    gpuErrchk( cudaEventCreate(&start) );
	    gpuErrchk( cudaEventCreate(&start1) );
	    gpuErrchk( cudaEventCreate(&stop) );
	    gpuErrchk( cudaEventCreate(&stop1) );
	    
	    // Record the start event
	    gpuErrchk( cudaEventRecord(start, NULL) );

	    // Launch Cuda Kernel
		dim3 block(BLOCK_SIZE, 1);
	    //dim3 grid(MAX_CONCURR_BL(BLOCK_SIZE), 1);
	    dim3 grid(1, 1);
	    
	    //GReset<<< grid, block >>>();
	    
	    //stConn<<< grid, block >>>(Dvertex, Dedges, Ddist_Col, N, nof_distNodes, DMatrix);
	    BFS_BlockKernel<<< grid, block, SM_BYTE_PER_BLOCK>>>(Dvertex, Dedges, Ddist_Col, Dsources, nof_distNodes);


		// Print matrix
		// gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMATRIX, cudaMemcpyDeviceToHost) );
		// for (int i = 0; i < nof_distNodes; ++i)
		// {
		// 	printf("| ");
		// 	for (int j = 0; j < nof_distNodes; ++j)
		// 		printf("%d ", matrix[nof_distNodes*i+j]);
		// 	printf("|\n");
		// }
		//printf("matrix completed\n");
		
	    gpuErrchk( cudaEventRecord(stop, NULL) );
	    gpuErrchk( cudaEventSynchronize(stop) );

	    bool connect = false;
	    if(!DEBUG){
		    gpuErrchk( cudaEventRecord(start1, NULL) );

		    // Copy result vector from device to host
			gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMATRIX, cudaMemcpyDeviceToHost) );
			connect = MatrixBFS(matrix, nof_distNodes, 0, 1, Queue);
			gpuErrchk( cudaEventRecord(stop1, NULL) );
		    gpuErrchk( cudaEventSynchronize(stop1) );
	    }


	    // Claculate elapsed time
	    float msecTotal = 0.0f;
	    float msecTotal1 = 0.0f;
	    gpuErrchk( cudaEventElapsedTime(&msecTotal, start, stop) );
	    gpuErrchk( cudaEventElapsedTime(&msecTotal1, start1, stop1) );

	    if(!DEBUG){
			printf("#%d:\tst-Connectivity from %d\t   to %d\tis %c[%d;%dm%s%c[%dm\t\tElapsed time = %c[%d;%dm%.1f%c[%dm ms\n", 
															test, source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
															27, 0, 27, 0, 31, msecTotal + msecTotal1, 27, 0);
			printf("-----------------------------------\n");	
	    }
		par_times[test] = msecTotal;
		seq_times[test] = msecTotal1;
	}
	double sum_par = 0;
	double sum_seq = 0;
	for (int i = 1; i < N_TEST; ++i){
		sum_par += par_times[i];
		sum_seq += seq_times[i];
	}
	if(!DEBUG){
		printf("\nN: %d\n", nof_distNodes);
		printf("AVG TIME \t\t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, (sum_par + sum_seq) / (N_TEST-1), 27, 0);
		printf("AVG PARALLEL TIME \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_par / (N_TEST-1), 27, 0);
		printf("AVG MATRIX BFS TIME \t: %c[%d;%dm%.1f%c[%dm ms\n\n", 27, 0, 31, sum_seq / (N_TEST-1), 27, 0);
	}

	//free memoria device
	cudaFree(Dvertex);
    cudaFree(Dedges);
    cudaFree(DMatrix);
    cudaFree(Ddist_Col);
    cudaFree(Ddistance);

    free(Queue);
	free(matrix);
	free(Dist_Col);
	free(Distance);
	free(vertex);
	free(edges);
	free(graph);

	return 0;
}