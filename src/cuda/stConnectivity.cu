#pragma once

#include "KernelStConnectivity_WEfficient.cu"
#include "stConn.h"

using namespace std;



void doSTCONN(Graph graph, int N, int E, int Nsources){
	
	/***    CALCULATE SIZES    ***/
	//int Nsources = atoi(argv[2]);

	size_t sizeE 	  = E * sizeof(int);
	size_t sizeN 	  = N * sizeof(int2);
	size_t sizeN1 	  = (N+1) * sizeof(int);
	size_t sizeN3 	  = Nsources * sizeof(int);
	size_t sizeMatrix = Nsources * Nsources * sizeof(bool);

	

	/***    ALLOCATE HOST MEMORY    ***/
	int2 *Distance 	= (int2*)calloc(N, sizeof(int2));	
	int *sources 	= (int*)calloc(Nsources, sizeof(int));
	int *Queue 		= (int*)calloc(Nsources, sizeof(int));	
	bool *matrix 	= (bool*)calloc(Nsources * Nsources, sizeof(bool));


	/***    ALLOCATE DEVICE MEMORY    ***/
	int *Dedges;
	int *Dvertex;
	int *Dsources;
	int2 *Ddistance;
	bool *DMatrix;

	gpuErrchk( cudaMalloc((void **) &Dvertex, 	sizeN1) );
	gpuErrchk( cudaMalloc((void **) &Dedges, 	sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, 	sizeMatrix) );
	gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN) );
	gpuErrchk( cudaMalloc((void **) &Dsources, 	sizeN3) );

	/***    PRINT CONFIG INFO    ***/
/*
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

	/***    SERVICE VARIABLES    ***/
 	vector<double> mean_times(3);
 	vector<double> par_times(N_TEST);
 	vector<double> seq_times(N_TEST);
 	vector<long double> Percentual(N_TEST);
 	int percCounter = 0;
 	int NULLCounter = 0;
 	long double perc = 0.0;


	srand (time(NULL));
 	for (int test = 0; test < N_TEST; test++)
	{

		/***    CHOOSE RANDOM SOURCE, DEST AND EXTRA-SOURCES    ***/
		int source = rand() % N;
		int target = rand() % N;

	    ChooseRandomNodes(sources, graph.OutNodes, N, Nsources, source, target);


	    /***    STRUCTURES INITIALIZATION    ***/
	    for (int i = 0; i < N; i++){
	    	Distance[i].x = INT_MAX;
	    	Distance[i].y = INT_MAX;
	    }

	    for (int i = 0; i < Nsources; i++){
			int j = sources[i];
			Distance[j].x = 0;
			Distance[j].y = i;
		}
		
		int VisitedNodes = 0;
	    
	    /***    MEMCOPY HOST_TO_DEVICE    ***/
	    gpuErrchk( cudaMemcpy(Dvertex, graph.OutNodes, sizeN1, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Dedges, graph.OutEdges, sizeE, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(DMatrix, matrix, sizeMatrix, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Ddistance, Distance, sizeN, cudaMemcpyHostToDevice) );
	    gpuErrchk( cudaMemcpy(Dsources, sources, sizeN3, cudaMemcpyHostToDevice) );		
		if(ATOMIC)
			gpuErrchk( cudaMemcpyToSymbol(GlobalCounter, &VisitedNodes, sizeof(int)) );
	    
		
		/***    ALLOCATE CUDA EVENT FOR TIMING    ***/
	    cudaEvent_t start;
	    cudaEvent_t stop;

	    gpuErrchk( cudaEventCreate(&start) );
	    gpuErrchk( cudaEventCreate(&stop) );
	    gpuErrchk( cudaEventRecord(start, NULL) );


		/***    LAUNCH KERNEL    ***/
		dim3 block(BLOCK_SIZE, 1);
	    dim3 grid(Nsources, 1);
	    BFS_BlockKernel<<< grid, block, SMem_Per_Block(BLOCK_SIZE)>>>(Dvertex, Dedges, Dsources, Ddistance, /*Dcolor,*/ DMatrix, Nsources);


	    /***    MEMCOPY DEVICE_TO_HOST    ***/
	    if(!BFS)
			gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );


	    /***    RECORD STOP TIME    ***/
	    gpuErrchk( cudaEventRecord(stop, NULL) );
	    gpuErrchk( cudaEventSynchronize(stop) );


	    /***    COPY EXITFLAG FROM DEVICE    ***/
	    if(!ATOMIC){
		    int Flag = 0;
		    gpuErrchk( cudaMemcpyFromSymbol(&Flag, exitFlag, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		    if(Flag){
		    	NULLCounter = N_TEST - test;
		    	break;
		    }
	    }


	    /***    CHECK VISIT PERCENTAGE IF IT FAILS    ***/
	    if(ATOMIC){
			gpuErrchk( cudaMemcpyFromSymbol(&VisitedNodes, GlobalCounter, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			VisitedNodes += Nsources;
			perc = ((long double)VisitedNodes / (long double)N) * 100;
			if(VisitedNodes < N){
				printf("---------------WARNING: BFS NOT COMPLETE---------------\t\t\t\t%.2Lf%\n", perc);
				Percentual[percCounter] = perc;
				percCounter++;
			}
		}

		/*gpuErrchk( cudaMemcpy(Distance, Ddistance, sizeN, cudaMemcpyDeviceToHost) );
		for (int i = 0; i < N; i++){
	    	if (Distance[i].x == INT_MAX)
	    		printf("---WARNING--- Nodo %d non visitato!!\n", i);
	    }*/


		

		/***    PRINT MATRIX    ***/
		// for (int i = 0; i < Nsources; ++i){
		// 	printf("| ");
		// 	for (int j = 0; j < Nsources; ++j)
		// 		printf("%d ", matrix[Nsources*i+j]);
		// 	printf("|\n");
		// }


	    float msecTotal = 0.0f;
	    float msecTotal1 = 0.0f;
	    bool connect = false;


	    /***    MATRIX VISIT ON HOST    ***/
	    if(!BFS){
		    Timer<HOST> TM;
		    TM.start();		    
			connect = MatrixBFS(matrix, Nsources, 0, 1, Queue);
			TM.stop();	    	
		    msecTotal1 = TM.duration();
	    }


		/***    CALCULATE ELAPSED TIME    ***/
	    gpuErrchk( cudaEventElapsedTime(&msecTotal, start, stop) );
	    
	    if(DEBUG)
			printf("#%d:\tsource: %d    \ttarget: %d   \tresult: %c[%d;%dm%s%c[%dm   \tElapsed time = %c[%d;%dm%.1f%c[%dm ms\n", 
															test, source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
															27, 0, 27, 0, 31, msecTotal + msecTotal1, 27, 0);
	    /*if( ATOMIC && perc < 100 ){
			par_times[test] = 0;
			seq_times[test] = 0;
			NULLCounter++;
	    }
	    else{*/
			par_times[test] = msecTotal;
			seq_times[test] = msecTotal1;
	   // }
	}
	

	/***    EVALUATE MEAN TIMES    ***/
	if(N_TEST > 1)
	{
		double sum_par = 0;
		double sum_seq = 0;
		for (int i = 1; i < N_TEST; ++i){
			sum_par += par_times[i];
			sum_seq += seq_times[i];
		}
		//printf("\nN: %d\n", Nsources);
		printf("# Completed task: %d on %d\n", N_TEST - NULLCounter, N_TEST);
		printf("AVG TIME \t\t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, (sum_par + sum_seq) / (N_TEST-NULLCounter), 27, 0);
		printf("AVG PARALLEL TIME \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_par / (N_TEST-NULLCounter), 27, 0);
		printf("AVG MATRIX BFS TIME \t: %c[%d;%dm%.1f%c[%dm ms\n\n", 27, 0, 31, sum_seq / (N_TEST-NULLCounter), 27, 0);
	}


	/***    EVALUATE MEAN PERCENTAGE    ***/
	if(ATOMIC){
		double sum = 0;
		for(int i = 0; i < percCounter; i++){
			sum += Percentual[i];
			//printf("sum = %f", sum);
		}
		printf("\n\nAVG Percentual \t\t: %.2f%\n", sum / percCounter);
		printf("MIN Percentual \t\t: %.2f%\n", min(Percentual, percCounter));
		printf("MAX Percentual \t\t: %.2f%\n", max(Percentual, percCounter));
	}

	
	/***    FREE MEMORY    ***/
	cudaFree(Dvertex);
    cudaFree(Dedges);
    cudaFree(DMatrix);
    cudaFree(Ddistance);

    free(Queue);
	free(matrix);
	free(Distance);
	//free(vertex);
	//free(edges);
	//free(graph);
}





int main(int argc, char *argv[]){

	if(argc < 2)
	{
		printf("\nUsage: ./stConnectivity 'input_file'\n\n");
		return -1;
	}


	/***    READ GRAPH FROM FILE    ***/
	int N, E, nof_lines;
 	GDirection GraphDirection;  	// scelta dell'utente oppure vuote (direzione di default estratta dal file) valori possibili = DIRECTED, UNDIRECTED
 	GraphDirection = UNDIRECTED;  			//DIRECTED = 0, UNDIRECTED = 1, UNDEFINED = 2

 	readGraph::readGraphHeader(argv[1], N, E, nof_lines, GraphDirection);
    Graph graph(N, E, GraphDirection);
    readGraph::readSTD(argv[1], graph, nof_lines);

    if(argc > 2)
    {
    	int Nsources = atoi(argv[2]);
    	printf("Launch stConnectivity with %d sources\n", Nsources);
		doSTCONN(graph, N, E, Nsources);
    }
    else
    {
		for (int i = 0; i < LENGTH; ++i)
		{
			printf("Launch stConnectivity with %d sources\n", parameters[i]);
			doSTCONN(graph, N, E, parameters[i]);
		}
    }

	

	return 0;
}
