//#pragma once

#include "Kernel_STCONN.cu"
#include "stConn.h"
#include "statistic.h"


void doSTCONN(Graph graph, int N, int E, double Treshold){
	
	/***    CALCULATE SIZES    ***/
	size_t sizeE 	  = E * sizeof(int);
	size_t sizeN 	  = N * sizeof(int2);
	size_t sizeN1 	  = (N+1) * sizeof(int);
	size_t sizeSrcs	  = 2 * sizeof(int);
	size_t sizeMatrix = MAX_SIZE * MAX_SIZE * sizeof(bool);
	size_t sizeBMask  = N * sizeof(bool);
	size_t sizeCurand = MAX_CONCURR_BL(BLOCK_SIZE) * sizeof(curandState);
	//size_t sizeBMask  = (N + (1024 * BLOCK_SIZE)) * sizeof(bool);
	

	/***    ALLOCATE HOST MEMORY    ***/
	int2 *Distance 	= (int2*)calloc(N, sizeof(int2));	
	int *sources 	= (int*)calloc(2, sizeof(int));
	int *Queue 		= (int*)calloc(MAX_SIZE, sizeof(int));	
	bool *matrix 	= (bool*)calloc(MAX_SIZE * MAX_SIZE, sizeof(bool));
	

	/***    ALLOCATE DEVICE MEMORY    ***/
	int *Dedges;
	int *Dvertex;
	int *Dsource;
	bool *DBitMask;
	int2 *Ddistance;
	bool *DMatrix;
	curandState* devStates;

    gpuErrchk( cudaMalloc((void **) &devStates, sizeCurand) );
	gpuErrchk( cudaMalloc((void **) &Dvertex, 	sizeN1) );
	gpuErrchk( cudaMalloc((void **) &Dedges, 	sizeE) );
	gpuErrchk( cudaMalloc((void **) &DMatrix, 	sizeMatrix) );
	gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN) );
	gpuErrchk( cudaMalloc((void **) &Dsource, 	sizeSrcs) );
	gpuErrchk( cudaMalloc((void **) &DBitMask, 	sizeBMask) );


	/***    SERVICE VARIABLES    ***/
	std::vector<double> mean_times(3);
	std::vector<double> par_times(N_TEST);
	std::vector<double> BOT_times(N_TEST);
	std::vector<double> seq_times(N_TEST);
	std::vector<long double> Percentual(N_TEST);
	int 	connectCnt = 0;
	int    percCounter = 0;
	int  unfinishedCnt = 0;
	long double   perc = 0.0;


	srand (time(NULL));
	for (int test = 0; test < N_TEST; ++test)
	{

		/***    CHOOSE RANDOM SOURCE, DEST AND EXTRA-SOURCES    ***/
		int VisitedEdges = 0;
		int source = rand() % N;
		int target = rand() % N;
		while(target == source)		target = rand() % N;


		/***    STRUCTURES INITIALIZATION    ***/
		for (int i = 0; i < N; ++i)
			Distance[i] = make_int2(INT_MAX, INT_MAX);

		sources[0] = source;
		sources[1] = target;
		Distance[source] = make_int2(0, 0);
		Distance[target] = make_int2(0, 1);


		/***    MEMCOPY HOST_TO_DEVICE    ***/
		gpuErrchk( cudaMemcpyToSymbol(GlobalCounter, &VisitedEdges, sizeof(int)) );
		gpuErrchk( cudaMemcpy(Dvertex, graph.OutNodes, sizeN1, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Dedges, graph.OutEdges, sizeE, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(DMatrix, matrix, sizeMatrix, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Ddistance, Distance, sizeN, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Dsource, sources, sizeSrcs, cudaMemcpyHostToDevice) );		
		gpuErrchk( cudaMemset(DBitMask, 0, sizeBMask) );

		
		/***    INITIALIZE TIMERS    ***/
		Timer<HOST> TM;
		Timer<DEVICE> TM_TD;
		Timer<DEVICE> TM_BU;
		float msecPAR = 0.0f;
		float msecSEQ = 0.0f;
		float msecBOT = 0.0f;
		bool  connect = false;


		/***    LAUNCH RESET KERNEL    ***/
		setup_curand<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>(devStates);

			int totSouce = 0;
			
		/***    LAUNCH STCONN TOP-DOWN KERNEL    ***/
		TM_TD.start();
		TopDown_Kernel<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
					(Dvertex, Dedges, Dsource, Ddistance, DMatrix, DBitMask, Treshold, E, N, devStates);
		TM_TD.stop();
		msecPAR = TM_TD.duration();

			cudaMemcpyFromSymbol(&totSouce, color, sizeof(int), 0, cudaMemcpyDeviceToHost);
			printf("Total sources = %d\n", totSouce);

		/***    CHECK VISIT PERCENTAGE    ***/
		gpuErrchk( cudaMemcpyFromSymbol(&VisitedEdges, GlobalCounter, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		perc = ((long double)VisitedEdges / (long double)E) * 100.0;
		

		if(BOTTOM_UP && perc > Treshold*100 )
		{

			//int FrontierSize1 = 0;
			//int level = 0;
			int FrontierSize = 1;
			int zero = 0;
			int result = 0;
			int BMsize = ceil((double)N / 4.0);

			TM_BU.start();
			while( FrontierSize )
			{
				//gpuErrchk( cudaMemcpyToSymbol(BottomUp_FrontSize1, &zero, sizeof(int), 0, cudaMemcpyHostToDevice) );
				gpuErrchk( cudaMemcpyToSymbol(BottomUp_FrontSize, &zero, sizeof(int), 0, cudaMemcpyHostToDevice) );
				BottomUp_Kernel<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
								(Dvertex, Dedges, DBitMask, BMsize, N);
				gpuErrchk( cudaMemcpyFromSymbol(&FrontierSize, BottomUp_FrontSize, sizeof(int), 0, cudaMemcpyDeviceToHost) );
				//gpuErrchk( cudaMemcpyFromSymbol(&FrontierSize1, BottomUp_FrontSize1, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			}
			TM_BU.stop();
			msecBOT = TM_BU.duration();



			cudaMemcpyToSymbol(VisitResult, &zero, sizeof(int), 0, cudaMemcpyHostToDevice);
			CheckVisit<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>(DBitMask, N);
			cudaMemcpyFromSymbol(&result, VisitResult, sizeof(int), 0, cudaMemcpyDeviceToHost);
			if(result != 0)
				printf("\t!!! There are %d nodes not visited !!!\n", result);
			
		}
		else
		{
			printf("---------------WARNING: TRESHOLD NOT REACHED---------------\t\t %d on %d\t%.2Lf%%\n", VisitedEdges, E, perc);
			Percentual[percCounter++] = perc;
		}


		/***    MATRIX STCONN ON HOST    ***/
		#if !BFS
			TM.start();		    
			
			//gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
			//connect = MatrixBFS(matrix, MAX_SIZE, 0, 1, Queue);

			TM.stop();	    	
			msecSEQ = TM.duration();
		#endif


		/***    PRINT RESULTS    ***/
		#if (!BFS)
			printf("#%d:\tsource: %d     \ttarget: %d      \tresult: %c[%d;%dm%s%c[%dm\t\ttime = %c[%d;%dm%.1f%c[%dm ms\n", 
															test+1, source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
															27, 0, 27, 0, 31, msecPAR + msecSEQ + msecBOT, 27, 0);
		#endif


		/***    SAVE TIMES FOR STATISTICS EVALUATION    ***/
		par_times[test] = msecPAR;
		seq_times[test] = msecSEQ;
		BOT_times[test] = msecBOT;

		if( perc < 100 )
			unfinishedCnt++;
	}
	

	/***    EVALUATE MEAN TIMES AND PERCENTAGE    ***/
	computeElapsedTime( par_times, seq_times, BOT_times, connectCnt);
	computeMeanPercentage(Percentual, percCounter);
	
	/***    FREE DEVICE MEMORY    ***/
	cudaFree(Ddistance);
	cudaFree(Dvertex);
	cudaFree(DMatrix);
	cudaFree(Dedges);
	cudaFree(DBitMask);
	cudaFree(Dsource);

	/***    FREE HOST MEMORY    ***/
	free(Distance);
	free(matrix);
	free(Queue);
	free(sources);
}



/*
* Read command line parameters
*/
void Parameters(int argc, char* argv[], GDirection &GDir, double& Treshold) {
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
        else if ( /*i + 1 < argc &&*/ parameter.compare("-n") == 0 && 
        		std::string(argv[i + 1]).find_first_not_of("0123456789.") == std::string::npos )
        {
            std::istringstream ss(argv[++i]);
            ss >> Treshold;
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
	double Treshold = 0.0;
	GDirection GraphDirection;		//DIRECTED = 0, UNDIRECTED = 1, UNDEFINED = 2
	Parameters(argc, argv, GraphDirection, Treshold);
	readGraph::readGraphHeader(argv[1], N, E, nof_lines, GraphDirection);
	Graph graph(N, E, GraphDirection);
	readGraph::readSTD(argv[1], graph, nof_lines);
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

	/***    LAUNCH ST-CONN FUNCTION    ***/
	if(Treshold != 0)
	{
		printf("\n----------Launch stConnectivity with treshold %.2f%%----------\n\n", Treshold*100);
		doSTCONN(graph, N, E, Treshold);
	}
	else
	{
		for (int i = 0; i < LENGTH; ++i)
		{
			printf("\n----------Launch stConnectivity with treshold %.2f%%----------\n\n", Treshold*100);
			doSTCONN(graph, N, E, Treshold);
		}
	}
	return 0;
}