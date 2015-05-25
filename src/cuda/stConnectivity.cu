#include <KernelStConnectivity.hpp>

using namespace std;


int main(int argc, char *argv[]){

	// If wrong number of argumants print usage
	if(argc < 3)
	{
		printf("\nUsage: ./stConnectivity 'input_file' '#_distinguished_nodes'\n\n");
		return -1;
	}

	// Read graph from file
	int x,y;
	edge *graph; 
	ifstream in (argv[1]);
	if(in.is_open())
	{
		//printf("Reading graph from file...\n");
		// Read graph parameters and allocate memory\n
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
	//printf("Graph sorted\n");

	// Memory allocation and initialization
	int *vertex = (int*)calloc(N+1, sizeof(int));
	int *edges = (int*)calloc(E, sizeof(int));
	//printf("Memory allocated\n");

	
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
	}

	// It is convenient to add a N-th element
	vertex[N] = E;

	//printf("CSR structure created\n");
	printf("-----------------------------------\n");

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
	
// printf("BLOCK_SIZE = %d\n", BLOCK_SIZE);
// int DD[13] = {10,50,100,500,1000,1500,2000,2500,3000,4000,5000,6000,7000};
// for (int s = 0; s < 13; ++s)
// {
// nof_distNodes = DD[s];

	// calculate size for allocation
	size_t sizeN = (N+1) * sizeof(int);
	size_t sizeE = E * sizeof(int);
	size_t sizeN2 = N * sizeof(int2);
	size_t sizeN3 = nof_distNodes * sizeof(int);
	size_t sizeMATRIX = nof_distNodes * nof_distNodes * sizeof(bool);

	// Allocate adj matrix in row major
	bool *matrix = (bool*)calloc(nof_distNodes * nof_distNodes, sizeof(bool));
	int2 *Dist_Col = (int2*)calloc(N, sizeof(int2));
	int *Distance = (int*)calloc(nof_distNodes, sizeof(int));	

	// Allocate device memory
	int *Dedges;
	int *Dvertex;
	bool *DMatrix;
	int2 *Ddist_Col;
	bool *DnewLevel;
	int *Ddistance;

	gpuErrchk( cudaMalloc((void **) &Dvertex, sizeN) );
	gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN3) );
	gpuErrchk( cudaMalloc((void **) &Dedges, sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, sizeMATRIX) );
	gpuErrchk( cudaMalloc((void **) &Ddist_Col, sizeN2) );
	// printf("Device memory allocated\n");
	// printf("-----------------------------------\n");

 
 	vector<double> par_times(N_TEST);
 	vector<double> seq_times(N_TEST);
 	for (int loop = 0; loop < 3; ++loop)
 	{
	 	for (int test = 0; test < N_TEST; test++)
		{
			srand (time(NULL));
			int source = rand() % N;
			int target = rand() % N;

			// choose nof_distNodes distinguished nodes with source and target in it and return a vector of them	
		    vector< int > id = ChooseNodes(vertex, N, nof_distNodes, source, target);
		    //printf("Distinguished nodes chosen\n");

		    // inizializzazione dei valori a INT_MAX
		    for (int i = 0; i < N; i++)
		    {
		    	Dist_Col[i].x = INT_MAX;
				Dist_Col[i].y = INT_MAX;
		    }

		    // inizializzazione dei valori dei nodi distinti a Distance 0 e Color id
		    for (int i = 0; i < nof_distNodes; i++)
			{
				//printf("choosen node : %d\n", id[i]);
				int j = id[i];
				Dist_Col[j].x = 0;
				Dist_Col[j].y = i;
				Distance[i] = INT_MAX;
			}
			Distance[0] = 0;
			bool NextLevel[4] = {0,0,0,0};
			
			//printf("Valori inizializzati\n");

			// for (int i = 0; i < 4; ++i)
			// {
			// 	NextLevel[i] = false;
			// }

		    // Copy host memory for vertex, edges and results vectors to device
		    gpuErrchk( cudaMemcpy(Dvertex, vertex, sizeN, cudaMemcpyHostToDevice) );
		    gpuErrchk( cudaMemcpy(Dedges, edges, sizeE, cudaMemcpyHostToDevice) );
		    gpuErrchk( cudaMemcpy(DMatrix, matrix, sizeMATRIX, cudaMemcpyHostToDevice) );
		    gpuErrchk( cudaMemcpy(Ddist_Col, Dist_Col, sizeN2, cudaMemcpyHostToDevice) );
		    gpuErrchk( cudaMemcpy(Ddistance, Distance, sizeN3, cudaMemcpyHostToDevice) );
		    gpuErrchk( cudaMemcpyToSymbolAsync(devNextLevel, NextLevel, sizeof(bool)*4, 0,cudaMemcpyHostToDevice) );
		    //printf("Memcopy executed\n");
			//printf("-----------------------------------\n");

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
		    dim3 grid(MAX_CONCURR_BL(BLOCK_SIZE), 1);
		    
		    GReset<<< grid, block >>>();
		    
		    stConn<<< grid, block >>>(Dvertex, Dedges, Ddist_Col, N, nof_distNodes, DMatrix);

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
		    // Wait for the stop event to complete
		    gpuErrchk( cudaEventSynchronize(stop) );
		    gpuErrchk( cudaEventRecord(start1, NULL) );

			bool connect = false;

		    // Copy result vector from device to host
			gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMATRIX, cudaMemcpyDeviceToHost) );
			connect = MatrixBFS(matrix, nof_distNodes, 0, 1);

			gpuErrchk( cudaEventRecord(stop1, NULL) );
		    // Wait for the stop event to complete
		    gpuErrchk( cudaEventSynchronize(stop1) );

		    // Claculate elapsed time
		    float msecTotal = 0.0f;
		    float msecTotal1 = 0.0f;
		    gpuErrchk( cudaEventElapsedTime(&msecTotal, start, stop) );
		    gpuErrchk( cudaEventElapsedTime(&msecTotal1, start1, stop1) );
		    //msecTotal /= 1000;
		    //msecTotal1 /= 1000;

			// printf("nodes = %d \tTime = %.3f s\n",nof_distNodes, msecTotal);
			// printf("Result for st-Connectivity from %d to %d is %c[%d;%dm%s%c[%dm\n", source, target, 27, 0, 
			// 												31 + connect,(connect ? "true" : "false"), 27, 0);
			// printf("Elapsed time \t= %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, msecTotal + msecTotal1, 27, 0);
			// // printf("Parallel time \t= %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, msecTotal, 27, 0);
			// // printf("MatrixBFS time \t= %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, msecTotal1, 27, 0);
			// printf("-----------------------------------\n");
			par_times[test] = msecTotal;
			seq_times[test] = msecTotal1;
		}
		double sum_par = 0;
		double sum_seq = 0;
		for (int i = 0; i < N_TEST; ++i){
			sum_par += par_times[i];
			sum_seq += seq_times[i];
		}
		printf("\nAVG TIME \t\t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, (sum_par + sum_seq) / N_TEST, 27, 0);
		printf("AVG PARALLEL TIME \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_par / N_TEST, 27, 0);
		printf("AVG MATRIX BFS TIME \t: %c[%d;%dm%.1f%c[%dm ms\n\n", 27, 0, 31, sum_seq / N_TEST, 27, 0);
	}
	//free memoria device
	cudaFree(Dvertex);
    cudaFree(Dedges);
    cudaFree(DMatrix);
    cudaFree(Ddist_Col);
    cudaFree(Ddistance);

	free(matrix);
	free(Dist_Col);
	free(Distance);
//}
	free(vertex);
	free(edges);
	free(graph);

	return 0;
}