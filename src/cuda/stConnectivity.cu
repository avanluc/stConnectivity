#pragma once

#include "Kernel_StConnectivity.cu"
#include "stConn.h"
#include "statistic.h"


void doSTCONN(Graph graph, int N, int E, int Nsources){
	
	/***    CALCULATE SIZES    ***/
	size_t sizeE 	  = E * sizeof(int);
	size_t sizeN 	  = N * sizeof(int2);
	size_t sizeN1 	  = (N+1) * sizeof(int);
	size_t sizeSrcs	  = Nsources * sizeof(int);
	size_t sizeMatrix = Nsources * Nsources * sizeof(bool);

	

	/***    ALLOCATE HOST MEMORY    ***/
	int2 *Distance 	= (int2*)calloc(N, sizeof(int2));	
	int *sources 	= (int*)calloc(Nsources, sizeof(int));
	int *Queue 		= (int*)calloc(Nsources, sizeof(int));	
	int *Visited    = (int*)calloc(Nsources, sizeof(int));
	bool *matrix 	= (bool*)calloc(Nsources * Nsources, sizeof(bool));


	/***    ALLOCATE DEVICE MEMORY    ***/
	int *Dedges;
	int *Dvertex;
	int *Dsources;
	int *Dvisited;
	int2 *Ddistance;
	bool *DMatrix;

	gpuErrchk( cudaMalloc((void **) &Dvertex, 	sizeN1) );
	gpuErrchk( cudaMalloc((void **) &Dedges, 	sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, 	sizeMatrix) );
	gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN) );
	gpuErrchk( cudaMalloc((void **) &Dsources, 	sizeSrcs) );
	gpuErrchk( cudaMalloc((void **) &Dvisited, 	sizeSrcs) );


	/***    SERVICE VARIABLES    ***/
 	std::vector<double> mean_times(3);
 	std::vector<double> par_times(N_TEST);
 	std::vector<double> seq_times(N_TEST);
 	std::vector<long double> Percentual(N_TEST);
 	int connectCnt = 0;
 	int percCounter = 0;
 	int unfinishedCnt = 0;
 	long double perc = 0.0;


	srand (time(NULL));
 	for (int test = 0; test < N_TEST; ++test)
	{

		/***    CHOOSE RANDOM SOURCE, DEST AND EXTRA-SOURCES    ***/
		int source = rand() % N;
		int target = rand() % N;
		while(target == source)		target = rand() % N;

	    ChooseRandomNodes(sources, N, Nsources, source, target);


	    /***    STRUCTURES INITIALIZATION    ***/
	    for (int i = 0; i < N; ++i){
	    	Distance[i].x = INT_MAX;
	    	Distance[i].y = INT_MAX;
	    }

	    for (int i = 0; i < Nsources; ++i){
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
	    gpuErrchk( cudaMemcpy(Dsources, sources, sizeSrcs, cudaMemcpyHostToDevice) );		
	    gpuErrchk( cudaMemcpy(Dvisited, Visited, sizeSrcs, cudaMemcpyHostToDevice) );		
		#if ATOMIC
			gpuErrchk( cudaMemcpyToSymbol(GlobalCounter, &VisitedNodes, sizeof(int)) );
	    #endif
		
		/***    ALLOCATE CUDA EVENT FOR TIMING    ***/
	    cudaEvent_t start;
	    cudaEvent_t stop;

	    gpuErrchk( cudaEventCreate(&start) );
	    gpuErrchk( cudaEventCreate(&stop) );
	    gpuErrchk( cudaEventRecord(start, NULL) );


		/***    LAUNCH KERNEL    ***/
		BFS_BlockKernel<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
	    			(Dvertex, Dedges, Dsources, Ddistance, DMatrix, Nsources);


	    /***    MEMCOPY DEVICE_TO_HOST    ***/
	    #if !BFS
			gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
		#endif


	    /***    RECORD STOP TIME    ***/
	    gpuErrchk( cudaEventRecord(stop, NULL) );
	    gpuErrchk( cudaEventSynchronize(stop) );


	    /***    COPY EXITFLAG FROM DEVICE    ***/
	    #if !ATOMIC
		    int Flag = 0;
		    gpuErrchk( cudaMemcpyFromSymbol(&Flag, exitFlag, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		    if(Flag){
		    	printf("Shared Memory Exceded!!!\n");
		    	unfinishedCnt = N_TEST - test;
		    	break;
		    }
	    #endif


	    /***    CHECK VISIT PERCENTAGE IF IT FAILS    ***/
	    #if(ATOMIC)
			gpuErrchk( cudaMemcpyFromSymbol(&VisitedNodes, GlobalCounter, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			VisitedNodes += Nsources;
			perc = ((long double)VisitedNodes / (long double)N) * 100;
			if(VisitedNodes < N){
				printf("---------------WARNING: BFS NOT COMPLETE---------------\t\t\t\t%.2Lf%\n", perc);
				Percentual[percCounter] = perc;
				percCounter++;
			}
		#endif

		/***    PRINT MATRIX    ***/
		//PrintMatrix<bool>(matrix, Nsources);


	    float msecTotal = 0.0f;
	    float msecTotal1 = 0.0f;
	    bool connect = false;


	    /***    MATRIX VISIT ON HOST    ***/
	    #if !BFS
		    Timer<HOST> TM;
		    TM.start();		    
			connect = MatrixBFS(matrix, Nsources, 0, 1, Queue);
			TM.stop();	    	
		    msecTotal1 = TM.duration();
	    #endif


		/***    CALCULATE ELAPSED TIME    ***/
	    gpuErrchk( cudaEventElapsedTime(&msecTotal, start, stop) );
	    
	    #if (!BFS && DEBUG)
			printf("#%d:\tsource: %d     \ttarget: %d      \tresult: %c[%d;%dm%s%c[%dm   \t\ttime = %c[%d;%dm%.1f%c[%dm ms\n", 
															test, source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
															27, 0, 27, 0, 31, msecTotal + msecTotal1, 27, 0);
		#endif
		par_times[test] = msecTotal;
		seq_times[test] = msecTotal1;
		
		if( !connect )
			connectCnt++;

	    if( ATOMIC && perc < 100 )
			unfinishedCnt++;
	}
	

	/***    EVALUATE MEAN TIMES    ***/
	#if(N_TEST > 1)
		computeElapsedTime( par_times, seq_times, connectCnt);
	#endif

	/***    EVALUATE MEAN PERCENTAGE    ***/
	#if ATOMIC
		computeMeanPercentage(Percentual, percCounter);
	#endif
	
	/***    FREE DEVICE MEMORY    ***/
    cudaFree(Ddistance);
	cudaFree(Dvertex);
    cudaFree(DMatrix);
    cudaFree(Dedges);

	/***    FREE HOST MEMORY    ***/
	free(Distance);
	free(matrix);
    free(Queue);
}



/*
* Read command line parameters
*/
void Parameters(int argc, char* argv[], GDirection &GDir, int& Nsources, int& all) {
    std::string errString(
    "Syntax Error:\n\n stConnectivity <graph_path> [ <graph_direction> ] [ -n <number_of_sources>] [-A]\n\n\
    <graph_direction>:\n\
                    -D      force directed graph\n\
                    -U      force undirected graph");

    if (argc < 2)
        error(errString)
    GDir = UNDEFINED;
    for (int i = 2; i < argc; ++i)
    {
        std::string parameter = argv[i];

        if 		(parameter.compare("-D") == 0)
            GDir = DIRECTED;
        else if (parameter.compare("-U") == 0)
            GDir = UNDIRECTED;
        else if (parameter.compare("-A") == 0)
        	all = 1;
        else if ( /*i + 1 < argc &&*/ parameter.compare("-n") == 0 && 
        		std::string(argv[i + 1]).find_first_not_of("0123456789") == std::string::npos )
        {
            std::istringstream ss(argv[++i]);
            ss >> Nsources;
        }
        else
            error(errString)
    }
}



/*
* Main function
*/
int main(int argc, char *argv[]){

	/***    READ GRAPH FROM FILE    ***/
	int N, E, nof_lines;
	int Nsources = 0, all = 0;
 	GDirection GraphDirection;		//DIRECTED = 0, UNDIRECTED = 1, UNDEFINED = 2
 	Parameters(argc, argv, GraphDirection, Nsources, all);
 	readGraph::readGraphHeader(argv[1], N, E, nof_lines, GraphDirection);
    Graph graph(N, E, GraphDirection);
    readGraph::readSTD(argv[1], graph, nof_lines);

    float avgDeg = (float) E / N;
    graph.DegreeAnalisys();


    /***    PRINT CONFIG INFO    ***/
	std::cout << "\n----------------------KERNEL INFO---------------------" 			<< std::endl
    	 << "            Block dimension : " <<  BLOCK_SIZE 							<< std::endl
    	 << "      Max concurrent blocks : " <<  MAX_CONCURR_BL(BLOCK_SIZE) 			<< std::endl
    	 << "       Shared Memory per SM : " <<  SMem_Per_SM 							<< std::endl
    	 << "    Shared Memory per block : " <<  SMem_Per_Block(BLOCK_SIZE) 			<< std::endl
    	 << "Int Shared Memory per block : " <<  IntSMem_Per_Block(BLOCK_SIZE) 			<< std::endl
    	 << "         Frontier dimension : " <<  FRONTIER_SIZE 							<< std::endl
    	 << "         Int frontier limit : " <<  BLOCK_FRONTIER_LIMIT 					<< std::endl
		 << "--------------------------------------------------------" 	   << std::endl << std::endl;

    /*if( graph.getMaxDegree() >= BLOCK_FRONTIER_LIMIT){
    	std::cout << std::endl << "Graph max degree greater than FRONTIER_LIMIT" << std::endl;
    	return 0;
    } */

	/***    LAUNCH ST-CONN FUNCTION    ***/
    if (all)
    {
		for (int i = 0; i < LENGTH; ++i)
		{
			printf("\n----------Launch stConnectivity with %d sources----------\n\n", SOURCES[i]);
			doSTCONN(graph, N, E, SOURCES[i]);
		}
    }
    else if(Nsources != 0)
    {
    	std::vector<double> prob = probability(N, Nsources, avgDeg);
    	for (int i = 0; i < prob.size(); ++i)
    		std::cout << "P(" << i << ") = " << std::right << std::setw(6) << prob[i] * 100 << "%" << std::endl;
    	printf("\nLaunch stConnectivity with %d sources\n\n", Nsources);
		doSTCONN(graph, N, E, Nsources);
    }
    else
    {
    	printf("Evaluating appropriate sources number\n");
    	Nsources = EvaluateSourcesNum(avgDeg, N);
    	printf("Launch stConnectivity with %d sources\n\n", Nsources);
		doSTCONN(graph, N, E, Nsources);
    }
	return 0;
}
