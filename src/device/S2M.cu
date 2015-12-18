#include "S2M_kernels.cu"
#include "S2M.h"
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

	bool *DBitMask1;
	int *DQueue;


	/***    RESULTS VECTORS    ***/
	std::vector<double> mean_times(3);
	std::vector<double> topDown_times(N_TEST);
	std::vector<double> bottomUp_times(N_TEST);
	std::vector<double> bfs_times(N_TEST);
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

		gpuErrchk( cudaMalloc((void **) &DBitMask1,   sizeBMask) );
		gpuErrchk( cudaMalloc((void **) &DQueue,      MAX_SIZE * sizeof(int)) );


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

		gpuErrchk( cudaMemset(DBitMask1, false, sizeBMask) );
		gpuErrchk( cudaMemset(DQueue, 0, MAX_SIZE * sizeof(int)) );

		
		/***    INITIALIZE TIMERS    ***/
		Timer<DEVICE> TM_TOP;
		Timer<DEVICE> TM_BOT;
		// Timer<HOST> TM_MATRIX1;
		// Timer<HOST> TM_MATRIX2;
		Timer<DEVICE> TM_MATRIX1;
		Timer<DEVICE> TM_MATRIX2;
		float msecTOP  = 0.0;
		float msecBOT  = 0.0;
		float msecBFS_TOP = 0.0;
		float msecBFS_BOT = 0.0;



		/***    LAUNCH RESET KERNEL    ***/
		clean<BLOCK_SIZE><<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>();



		/***    LAUNCH STCONN TOP-DOWN KERNEL    ***/
		TM_TOP.start();
		TopDown_Kernel<BLOCK_SIZE, V_WARP_SIZE><<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
					(Dvertex, Dedges, Dsource, Ddistance, DMatrix, DBitMask, devStates, Treshold, E, N);
		gpuErrchk( cudaMemcpyFromSymbol(&tempSources, color, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		TM_TOP.stop();


		/***    CHECK VISIT PERCENTAGE AND SOURCES    ***/
		gpuErrchk( cudaMemcpyFromSymbol(&VisitedEdges, GlobalCounter, sizeof(int), 0, cudaMemcpyDeviceToHost) );
		perc = ((long double)VisitedEdges / (long double)E) * 100.0;
		totSrc += (tempSources-1);



		/***    MATRIX STCONN ON HOST    ***/
/*
		TM_MATRIX1.start();
		gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
		connect1 = MatrixBFS(matrix, tempSources, 0, 1, Queue);
		TM_MATRIX1.stop();
*/

		TM_MATRIX1.start();
		matrixVisit<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
					(DMatrix, DBitMask1, DQueue, totSrc, 0, 1);
		TM_MATRIX1.stop();
		gpuErrchk( cudaMemcpyFromSymbol(&connect1, matrixResult, sizeof(bool), 0, cudaMemcpyDeviceToHost) );



		msecTOP = TM_TOP.duration();
		msecBFS_TOP = TM_MATRIX1.duration();


		if(BOTTOM_UP && !connect1 && perc >= Treshold*100 )
		{
			int result = 0;
			int result1 = 0;
			int FrontierSize = 1;


			/***    LAUNCH STCONN BOTTOM_UP KERNEL    ***/
			TM_BOT.start();
			while( FrontierSize )
			{
				gpuErrchk( cudaMemcpyToSymbol(BottomUp_FrontSize, &ZERO, sizeof(int), 0, cudaMemcpyHostToDevice) );
				BottomUp_Kernel<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
								(Dvertex, Dedges, Ddistance, DMatrix, DBitMask, N);
				gpuErrchk( cudaMemcpyFromSymbol(&FrontierSize, BottomUp_FrontSize, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			}
			TM_BOT.stop();



			/***    CHECK IF THE VISIT IS COMPLETE    ***/
			CheckVisit<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>(DBitMask, Ddistance, N);

			gpuErrchk( cudaMemcpyFromSymbol(&result, VisitResult, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			gpuErrchk( cudaMemcpyFromSymbol(&result1, VisitResult1, sizeof(int), 0, cudaMemcpyDeviceToHost) );
			if(result1)		printf("\t!!! There are %d nodes without color !!!\n", result1 );
			if(result)		printf("\t!!! There are %d nodes not visited !!!\n", result);


			/***    MATRIX STCONN ON HOST    ***/
/*
			TM_MATRIX2.start();
			gpuErrchk( cudaMemcpy(matrix, DMatrix, sizeMatrix, cudaMemcpyDeviceToHost) );
			connect = MatrixBFS(matrix, tempSources, 0, 1, Queue);
			TM_MATRIX2.stop();
*/

			gpuErrchk( cudaMemset(DBitMask1, false, sizeBMask) );
			TM_MATRIX2.start();
			matrixVisit<<< MAX_CONCURR_BL(BLOCK_SIZE), BLOCK_SIZE, SMem_Per_Block(BLOCK_SIZE)>>>\
					(DMatrix, DBitMask1, DQueue, totSrc, 0, 1);
			TM_MATRIX2.stop();
			gpuErrchk( cudaMemcpyFromSymbol(&connect, matrixResult, sizeof(bool), 0, cudaMemcpyDeviceToHost) );

			

			msecBOT = TM_BOT.duration();
			msecBFS_BOT = TM_MATRIX2.duration();

			PrintResults(test, source, target, connect, msecTOP + msecBFS_TOP + msecBOT + msecBFS_BOT);
		}
		else if( perc < Treshold*100 )
		{
			printf("---------------WARNING: TRESHOLD NOT REACHED---------------\t\t\t%.2Lf%%\n", perc);
		}
		else
		{
			PrintResults(test, source, target, connect1, msecTOP + msecBFS_TOP);
		}



		#if SEQ_CHECK
			bool conn = stConnectivity(graph.OutNodes, graph.OutEdges, N, source, target);
			PrintResults(test, source, target, conn, -1.0);
		#endif



		/***    SAVE TIMES FOR STATISTICS EVALUATION    ***/
		topDown_times[test] = msecTOP;
		bottomUp_times[test] = msecBOT;
		bfs_times[test] = msecBFS_BOT + msecBFS_TOP;


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
	computeElapsedTime( topDown_times, bfs_times, bottomUp_times);
	printf("AVG SOURCES        \t: %d\n", totSrc/N_TEST);
	

}


