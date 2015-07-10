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
	printf("Transfert graph to CSR structure\n");
	
	// int source = atoi(argv[2]);
	// int target = atoi(argv[3]);
	// int Nsources = atoi(argv[4]);

	int Nsources = atoi(argv[2]);

	// calculate size for allocation
	size_t sizeE 	= E * sizeof(int);
	//size_t sizeN 	= N * sizeof(int);
	size_t sizeN 	= N * sizeof(int2);
	size_t sizeN1 	= (N+1) * sizeof(int);
	size_t sizeN3 	= Nsources * sizeof(int);
	size_t sizeMatrix = Nsources * Nsources * sizeof(bool);

	// Allocate adj matrix in row major
	// int *Color 		= (int*)calloc(N, sizeof(int));
	// int *Distance 	= (int*)calloc(N, sizeof(int));	
	int2 *Distance 	= (int2*)calloc(N, sizeof(int2));	
	int *sources 	= (int*)calloc(Nsources, sizeof(int));
	int *Queue 		= (int*)calloc(Nsources, sizeof(int));	
	bool *matrix 	= (bool*)calloc(Nsources * Nsources, sizeof(bool));

	// Allocate device memory
	int *Dedges;
	int *Dvertex;
	int *Dsources;
	int2 *Ddistance;
	// int *Ddistance;
	// int *Dcolor;
	bool *DMatrix;

	gpuErrchk( cudaMalloc((void **) &Dvertex, 	sizeN1) );
	gpuErrchk( cudaMalloc((void **) &Dedges, 	sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, 	sizeMatrix) );
	gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN) );
	//gpuErrchk( cudaMalloc((void **) &Dcolor, 	sizeN) );
	gpuErrchk( cudaMalloc((void **) &Dsources, 	sizeN3) );
/*
	printf("Device memory allocated\n");

	cout << "\n---------------------KERNEL INFO---------------------" 					<< endl
    	 << "            Block dimension : " <<  BLOCK_SIZE 							<< endl
    	 << "      Max concurrent blocks : " <<  MAX_CONCURR_BL(BLOCK_SIZE) 			<< endl
    	 << "   Number of current blocks : " <<  Nsources 								<< endl
    	 << "       Shared Memory per SM : " <<  SMem_Per_SM 							<< endl
    	 << "    Shared Memory per block : " <<  SMem_Per_Block(BLOCK_SIZE) 			<< endl
    	 << "Int Shared Memory per block : " <<  IntSMem_Per_Block(BLOCK_SIZE) 			<< endl
    	 << "                  F1_OFFSET : " <<  F1_OFFSET 								<< endl
    	 << "               F1 dimension : " <<  SMem_Per_Block(BLOCK_SIZE)-F1_OFFSET	<< endl
    	 << "         Frontier dimension : " <<  FRONTIER_SIZE 							<< endl
    	 << "       Block frontier limit : " <<  BLOCK_FRONTIER_LIMIT 					<< endl
		 << "-------------------------------------------------------" 			<< endl << endl;
*/
 	//vector< ONodes > ordered = OrderNodes(vertex, N);

 	vector<double> mean_times(3);
 	vector<double> par_times(N_TEST);
 	vector<double> seq_times(N_TEST);
 	vector<long double> Percentual(N_TEST);
 	int percCounter = 0;


	srand (time(NULL));
 	for (int test = 0; test < N_TEST; test++)
	{

		int source = rand() % N;
		int target = rand() % N;
		//printf("#%d:\tsource: %d    \ttarget: %d\n", test,source, target);

		// choose Nsources distinguished nodes with source and target in it and return a vector of them	
	    //ChooseNodes(sources, ordered, Nsources, source, target);
	    ChooseRandomNodes(sources, vertex, N, Nsources, source, target);

	    // inizializzazione dei valori a INT_MAX
	    for (int i = 0; i < N; i++){
	    	Distance[i].x = INT_MAX;
	    	Distance[i].y = INT_MAX;
	  		// Distance[i] = INT_MAX;
			// Color[i] = INT_MAX;
	    }

	    // inizializzazione dei valori dei nodi distinti a Distance 0 e Color id
	    for (int i = 0; i < Nsources; i++){
			int j = sources[i];
			Distance[j].x = 0;
			Distance[j].y = i;
			// Distance[j] = 0;
			// Color[j] = i;
			//printf("sources[%d] = %d\n",i,j );
		}
		
		int VisitedNodes = 0;
	    // Copy host memory for vertex, edges and results vectors to device
	    gpuErrchk( cudaMemcpy(Dvertex, vertex, sizeN1, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(DMatrix, matrix, sizeMatrix, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Dedges, edges, sizeE, cudaMemcpyHostToDevice) );
	    //gpuErrchk( cudaMemcpy(Dcolor, Color, sizeN, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Ddistance, Distance, sizeN, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Dsources, sources, sizeN3, cudaMemcpyHostToDevice) );
		
		gpuErrchk( cudaMemcpyToSymbol(GlobalCounter, &VisitedNodes, sizeof(int)) );
	    
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

		// if(SINGLE_BLOCK){
		//     dim3 grid(1, 1);
		//     BFS_BlockKernel<<< grid, block, SMem_Per_SM>>>(Dvertex, Dedges, Dsources, Ddistance, /*Dcolor,*/ DMatrix, Nsources);
		// }
		// else{
		    dim3 grid(Nsources, 1);
		    BFS_BlockKernel<<< grid, block, SMem_Per_Block(BLOCK_SIZE)>>>(Dvertex, Dedges, Dsources, Ddistance, /*Dcolor,*/ DMatrix, Nsources);
		// }

		// gpuErrchk( cudaMemcpyFromSymbol(&VisitedNodes, GlobalCounter, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		// VisitedNodes += Nsources;
		// long double perc = ((long double)VisitedNodes / (long double)N) * 100;
		// //cout << "              Visited Nodes : " << VisitedNodes << endl 
		// //	 << "                Graph Nodes : "<< N << endl;
		// //printf("  Graph Visitage Percentual : %.2Lf%\n", perc);
		// if(VisitedNodes < N){
		// 	printf("---------------WARNING: BFS NOT COMPLETE---------------\t\t\t\t%.2Lf%\n", perc);
		// 	Percentual[percCounter] = perc;
		// 	percCounter++;
		// }
		// else{
		// 	printf("-------------------------------------------------------\n");
		// }

		/*gpuErrchk( cudaMemcpy(Distance, Ddistance, sizeN, cudaMemcpyDeviceToHost) );
		for (int i = 0; i < N; i++){
	    	if (Distance[i].x == INT_MAX)
	    		printf("---WARNING--- Nodo %d non visitato!!\n", i);
	    }*/

		gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
		// for (int i = 0; i < Nsources; ++i)
		// {
		// 	printf("| ");
		// 	for (int j = 0; j < Nsources; ++j)
		// 		printf("%d ", matrix[Nsources*i+j]);
		// 	printf("|\n");
		// }
		//printf("matrix completed\n");
		
	    gpuErrchk( cudaEventRecord(stop, NULL) );
	    gpuErrchk( cudaEventSynchronize(stop) );
		gpuErrchk( cudaEventRecord(start1, NULL) );
	    
	    bool connect = false;
		connect = MatrixBFS(matrix, Nsources, 0, 1, Queue);

		gpuErrchk( cudaEventRecord(stop1, NULL) );
	    gpuErrchk( cudaEventSynchronize(stop1) );

	    // Claculate elapsed time
	    float msecTotal = 0.0f;
	    float msecTotal1 = 0.0f;
	    gpuErrchk( cudaEventElapsedTime(&msecTotal, start, stop) );
	    gpuErrchk( cudaEventElapsedTime(&msecTotal1, start1, stop1) );

	    //if(!DEBUG){
			// printf("#%d:\tsource: %d    \ttarget: %d   \tresult: %c[%d;%dm%s%c[%dm   \tElapsed time = %c[%d;%dm%.1f%c[%dm ms\n", 
			// 												test, source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
			// 												27, 0, 27, 0, 31, msecTotal + msecTotal1, 27, 0);
			//printf("Parallel Time : %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, msecTotal, 27, 0);

			//printf("-----------------------------------\n");	
	    //}
		par_times[test] = msecTotal;
		seq_times[test] = msecTotal1;
	}
	if(N_TEST > 1)
	{
		double sum_par = 0;
		double sum_seq = 0;
		for (int i = 1; i < N_TEST; ++i){
			sum_par += par_times[i];
			sum_seq += seq_times[i];
		}
		printf("\nN: %d\n", Nsources);
		printf("AVG TIME \t\t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, (sum_par + sum_seq) / (N_TEST-1), 27, 0);
		printf("AVG PARALLEL TIME \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_par / (N_TEST-1), 27, 0);
		printf("AVG MATRIX BFS TIME \t: %c[%d;%dm%.1f%c[%dm ms\n\n", 27, 0, 31, sum_seq / (N_TEST-1), 27, 0);
	}

	// double sum = 0;
	// for(int i = 0; i < percCounter; i++){
	// 	sum += Percentual[i];
	// 	printf("sum = %f", sum);
	// }
	// printf("\n\nAVG Percentual \t\t: %.2f%\n", sum / percCounter);
	// printf("MIN Percentual \t\t: %.2f%\n", min(Percentual, percCounter));
	// printf("MAX Percentual \t\t: %.2f%\n", max(Percentual, percCounter));

	//free memory
	cudaFree(Dvertex);
    cudaFree(Dedges);
    cudaFree(DMatrix);
    //cudaFree(Dcolor);
    cudaFree(Ddistance);

    free(Queue);
	free(matrix);
	//free(Color);
	free(Distance);
	free(vertex);
	free(edges);
	free(graph);

	return 0;
}