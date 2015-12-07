//#pragma once

#include "STCONN_kernels.cu"
#include "stConn.h"
#include "statistic.h"



void doSTCONN(Graph graph, int N, int E, double Treshold){
	
	/***    CALCULATE SIZES    ***/
	size_t sizeE      = E * sizeof(int);
	size_t sizeN      = N * sizeof(int2);
	size_t sizeSrcs	  = 2 * sizeof(int);
	size_t sizeBMask  = N * sizeof(bool);
	size_t sizeN1     = (N+1) * sizeof(int);
	size_t sizeMatrix = MAX_SIZE * MAX_SIZE * sizeof(bool);
	size_t sizeCurand = MAX_CONCURR_BL(BLOCK_SIZE) * sizeof(curandState);
	

	

	/***    DEVICE MEMORY POINTERS    ***/
	int *Dedges;
	int *Dvertex;
	int *Dsource;
	bool *DBitMask;
	int2 *Ddistance;
	bool *DMatrix;
	curandState* devStates;


	/***    RESULTS VECTORS    ***/
	std::vector<double> mean_times(3);
	std::vector<double> topDown_times(N_TEST);
	std::vector<double> bottomUp_times(N_TEST);
	std::vector<double> seq_times(N_TEST);
	std::vector<long double> Percentual(N_TEST);
	int totSrc = 0;

    
    /***    INITIALIZE CURAND SEEDS    ***/
    gpuErrchk( cudaMalloc((void **) &devStates, sizeCurand) );
	setup_curand<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>(devStates);


	srand (time(NULL));
	for (int test = 0; test < N_TEST; ++test)
	{
		/***    ALLOCATE HOST MEMORY    ***/
		int2 *Distance = (int2*)calloc(N, sizeof(int2));	
		int *sources   = (int*)calloc(2, sizeof(int));
		int *Queue     = (int*)calloc(MAX_SIZE, sizeof(int));	
		bool *matrix   = (bool*)calloc(MAX_SIZE * MAX_SIZE, sizeof(bool));


		/***    SERVICE VARIABLES    ***/
		int  VisitedEdges = 0;
		int   tempSources = 0;
		bool      connect = false;
		bool     connect1 = false;
		long double  perc = 0.0;


		/***    ALLOCATE DEVICE MEMORY    ***/
		gpuErrchk( cudaMalloc((void **) &Dvertex,   sizeN1) );
		gpuErrchk( cudaMalloc((void **) &Dedges,    sizeE) );
		gpuErrchk( cudaMalloc((void **) &DMatrix,   sizeMatrix) );
		gpuErrchk( cudaMalloc((void **) &Ddistance, sizeN) );
		gpuErrchk( cudaMalloc((void **) &Dsource,   sizeSrcs) );
		gpuErrchk( cudaMalloc((void **) &DBitMask,  sizeBMask) );


		/***    CHOOSE RANDOM SOURCE AND DEST    ***/
		int source = rand() % N;
		int target = rand() % N;
		while(target == source) 	target = rand() % N;


		/***    STRUCTURES INITIALIZATION    ***/
		for (int i = 0; i < N; ++i)
			Distance[i] = make_int2(INT_MAX, INT_MAX);

		sources[0] = source;
		sources[1] = target;
		Distance[source] = make_int2(0, 0);
		Distance[target] = make_int2(0, 1);


		/***    MEMCOPY HOST_TO_DEVICE    ***/
		gpuErrchk( cudaMemcpyToSymbol(VisitResult,   &ZERO, sizeof(int), 0, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpyToSymbol(VisitResult1,  &ZERO, sizeof(int), 0, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpyToSymbol(GlobalCounter, &ZERO, sizeof(int), 0, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpyToSymbol(color, &INIT_COLOR, sizeof(int), 0, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Dvertex, graph.OutNodes, sizeN1, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Dedges,  graph.OutEdges, sizeE,  cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Ddistance, Distance, sizeN, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(Dsource, sources, sizeSrcs, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemset(DMatrix, false, sizeMatrix) );
		gpuErrchk( cudaMemset(DBitMask, false, sizeBMask) );

		
		/***    INITIALIZE TIMERS    ***/
		Timer<HOST> TM_SEQ;
		Timer<HOST> TM_SEQ1;
		Timer<DEVICE> TM_TD;
		Timer<DEVICE> TM_BU;
		float msecTOP  = 0.0f;
		float msecBOT  = 0.0f;
		float msecSEQ_BOT  = 0.0f;
		float msecSEQ_TOP = 0.0f;



		/***    LAUNCH RESET KERNEL    ***/
		clean<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>();
	


		/***    LAUNCH STCONN TOP-DOWN KERNEL    ***/
		TM_TD.start();
		TopDown_Kernel<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
					(Dvertex, Dedges, Dsource, Ddistance, DMatrix, DBitMask, Treshold, E, N, devStates);
		TM_TD.stop();



		/***    MATRIX STCONN ON HOST    ***/
		TM_SEQ1.start();
		gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
		connect1 = MatrixBFS(matrix, MAX_SIZE, 0, 1, Queue);
		TM_SEQ1.stop();



		msecTOP = TM_TD.duration();
		msecSEQ_TOP = TM_SEQ1.duration();



		/***    CHECK VISIT PERCENTAGE AND SOURCES    ***/
		gpuErrchk( cudaMemcpyFromSymbol(&tempSources, color, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		gpuErrchk( cudaMemcpyFromSymbol(&VisitedEdges, GlobalCounter, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		perc = ((long double)VisitedEdges / (long double)E) * 100.0;
		totSrc += (tempSources-1);
		

		if(BOTTOM_UP && !connect1 && perc >= Treshold*100 )
		{
			int result = 0;
			int result1 = 0;
			int FrontierSize = 1;


			/***    LAUNCH STCONN BOTTOM_UP KERNEL    ***/
			TM_BU.start();
			while( FrontierSize )
			{
				gpuErrchk( cudaMemcpyToSymbol(BottomUp_FrontSize, &ZERO, sizeof(int), 0, cudaMemcpyHostToDevice) );
				BottomUp_Kernel<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
								(Dvertex, Dedges, Ddistance, DMatrix, DBitMask, N);
				gpuErrchk( cudaMemcpyFromSymbol(&FrontierSize, BottomUp_FrontSize, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			}
			TM_BU.stop();



			/***    CHECK IF THE VISIT IS COMPLETE    ***/
			CheckVisit<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>(DBitMask, Ddistance, N);

			gpuErrchk( cudaMemcpyFromSymbol(&result, VisitResult, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			gpuErrchk( cudaMemcpyFromSymbol(&result1, VisitResult1, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			if(result1)		printf("\t!!! There are %d nodes without color !!!\n", result1 );
			if(result)		printf("\t!!! There are %d nodes not visited !!!\n", result);


			/***    MATRIX STCONN ON HOST    ***/
			TM_SEQ.start();
			gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
			connect = MatrixBFS(matrix, MAX_SIZE, 0, 1, Queue);
			TM_SEQ.stop();
			

			msecBOT = TM_BU.duration();
			msecSEQ_BOT = TM_SEQ.duration();


			PrintResults(test, source, target, connect, msecTOP + msecSEQ_TOP + msecBOT + msecSEQ_BOT);
		}
		else if( perc < Treshold*100 )
		{
			printf("---------------WARNING: TRESHOLD NOT REACHED---------------\t\t\t%.2Lf%%\n", perc);
		}
		else
		{
			PrintResults(test, source, target, connect1, msecTOP + msecSEQ_TOP);
		}



		#if SEQ_CHECK
			bool conn = stConnectivity(graph.OutNodes, graph.OutEdges, N, source, target);
			PrintResults(test, source, target, conn, -1.0);
		#endif



		/***    SAVE TIMES FOR STATISTICS EVALUATION    ***/
		topDown_times[test] = msecTOP;
		bottomUp_times[test] = msecBOT;
		seq_times[test] = msecSEQ_BOT + msecSEQ_TOP;


		/***    FREE DEVICE MEMORY    ***/
		cudaFree(Ddistance);
		cudaFree(DBitMask);
		cudaFree(Dsource);
		cudaFree(DMatrix);
		cudaFree(Dvertex);
		cudaFree(Dedges);
		

		/***    FREE HOST MEMORY    ***/
		free(Distance);
		free(matrix);
		free(Queue);
		free(sources);
	}
	cudaFree(devStates);


	/***    EVALUATE MEAN TIMES AND PERCENTAGE    ***/
	computeElapsedTime( topDown_times, seq_times, bottomUp_times);
	printf("AVG SOURCES        \t: %d\n", totSrc/N_TEST);
	

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
			printf("\n----------Launch stConnectivity with treshold %.2f%%----------\n\n", TRESH[i]*100);
			doSTCONN(graph, N, E, TRESH[i]);
		}
	}
	return 0;
}