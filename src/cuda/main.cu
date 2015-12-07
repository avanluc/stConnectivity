#include "stConnectivity.cu"


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