#include "graph.h"

Graph::Graph(const int _V, const int _E, const GDirection GraphDirection) :
				V(_V), E(_E), Direction(GraphDirection) {
	try {
		Queue = new int[V];
		Visited.resize(V);

		OutNodes = new int[ V + 1 ];
		OutNodes2 = new int2[ V ];
		OutEdges = new int[ E ];
		OutDegree = new int[ V ]();
		COO_Edges = new int2[ E ];
		if (GraphDirection == UNDIRECTED)
			return;
		InNodes = new int[ V + 1 ];
		InEdges = new int[ E ];
		InDegree = new int[ V ]();
	}
	catch(std::bad_alloc& exc) {
  		error("OUT OF MEMORY: Graph too Large (Out Graph)!!");
	}
}

void Graph::ToCSR() {
	std::cout << "        COO To CSR...\t\t" << std::flush;

	for (int i = 0; i < COOSize; i++) {
		int index1 = COO_Edges[i].x;
		int index2 = COO_Edges[i].y;
		OutDegree[index1]++;
		if (Direction == UNDIRECTED)
			OutDegree[index2]++;
		else if (Direction == DIRECTED)
			InDegree[index2]++;
	}
	OutNodes[0] = 0;
	std::partial_sum(OutDegree, OutDegree + V, OutNodes + 1);
	for (int i = 0; i < V; i++)
		OutNodes2[i] = make_int2(OutNodes[i], OutNodes[i + 1]);

	int* TMP = new int[V]();
	for (int i = 0; i < COOSize; ++i) {
		int index1 = COO_Edges[i].x;
		int index2 = COO_Edges[i].y;
		OutEdges[ OutNodes[index1] + TMP[index1]++ ] = index2;
		if (Direction == UNDIRECTED)
			OutEdges[ OutNodes[index2] + TMP[index2]++ ] = index1;
	}

	if (Direction == DIRECTED) {
		InNodes[0] = 0;
		std::partial_sum(InDegree, InDegree + V, InNodes + 1);

		std::fill(TMP, TMP + V, 0);
		for (int i = 0; i < COOSize; ++i) {
			int index2 = COO_Edges[i].y;
			InEdges[ InNodes[index2] + TMP[index2]++ ] = COO_Edges[i].x;
		}
	}
	delete TMP;
	std::cout << "Complete!\n\n" << std::flush;
}

void  Graph::Dimacs10ToCOO() {
	std::cout << " Dimacs10th to COO...\t\t" << std::flush;
	int count_Edges = 0;
	for (int i = 0; i < V; i++) {
		for (int j = OutNodes[i]; j < OutNodes[i + 1]; j++) {
			int dest = OutEdges[j];
			bool flag = true;
			for (int t = OutNodes[dest]; flag && t < OutNodes[dest + 1]; t++) {
				if (t == i)
					flag = false;
			}
			if (flag)
				COO_Edges[count_Edges++] = make_int2(i, dest);
		}
	}
	COOSize = count_Edges;
	std::cout << "Complete!\n\n" << std::flush;
}

void Graph::print() {
	printExt::printArray(COO_Edges, COOSize, "COO Edges\n");
	printExt::printArray(OutNodes, V + 1, "OutNodes\t");
	printExt::printArray(OutEdges, E, "OutEdges\t");
	printExt::printArray(OutDegree, V, "OutDegree\t");
	if (Direction == UNDIRECTED)
		return;
	printExt::printArray(InNodes, V + 1, "InNodes\t\t");
	printExt::printArray(InEdges, E, "InEdges\t\t");
	printExt::printArray(InDegree, V, "InDegree\t");
}

void Graph::BfsInit(int source, int* _Distance) {
	left = 0, right = 1;

	std::fill(Visited.begin(), Visited.end(), false);
	Visited[source] = true;

	Distance = _Distance;
	std::fill(Distance, Distance + V, std::numeric_limits<int>::max());
	Distance[source] = 0;

	Queue[0] = source;
}

void Graph::bfs() {
	while (left < right) {
		int qNode = Queue[left++];

		for (int i = OutNodes[qNode]; i < OutNodes[qNode + 1]; ++i) {
			int dest = OutEdges[i];

			if (!Visited[dest]) {
				Visited[dest] = true;
				Distance[dest] = Distance[qNode] + 1;
				Queue[right++] = dest;
			}
		}
	}
}

int Graph::visitedNodes() {
	return right;
}

int Graph::visitedEdges() {
	if (right == V)
		return E;
	int sum = 0;
	for (int i = 0; i < right; ++i)
		sum += OutDegree[ Queue[i] ];
	return sum;
}

int Graph::getMaxDistance() {
	return Distance[Queue[right - 1]];
}

int Graph::getMaxDegree(){
	std::pair<int*,int*> minmax = std::minmax_element(OutDegree, OutDegree + V);
	return *minmax.second;
}

void Graph::bfsFrontier(std::vector<int>& Frontiers) {
	int oldDistance = 0;
	Frontiers.push_back(1);
	int edges = 0;
	int level = 0;

	while (left < right) {
		int qNode = Queue[left++];

		if (Distance[qNode] > oldDistance) {
			Frontiers.push_back(right - left + 1);
			oldDistance = Distance[qNode];
			//std::cout << std::endl << std::endl << "LEVEL " << oldDistance + 1 << std::endl;
			std::cout << std::endl << std::endl << "LEVEL " << oldDistance + 1 << " Edges " << edges << std::endl;
			edges = 0;
		}
		edges += OutDegree[qNode];
		//std::cout <<  qNode << ' ';

		for (int i = OutNodes[qNode]; i < OutNodes[qNode + 1]; ++i) {
			int dest = OutEdges[i];

			std::cout << dest << " ";
			if (!Visited[dest]) {
				Visited[dest] = true;
				Distance[dest] = Distance[qNode] + 1;
				Queue[right++] = dest;
			}
		}
	}
}

void Graph::DegreeAnalisys() {
	std::cout.imbue(std::locale(std::locale(), new fUtil::myseps));
	float avg = (float) E / V;
	float stdDev = fUtil::stdDeviation(OutDegree, V, avg);
	std::pair<int*,int*> minmax = std::minmax_element(OutDegree, OutDegree + V);
	//int zeroDegree = std::count (OutDegree, OutDegree + V, 0);
	//int oneDegree = std::count (OutDegree, OutDegree + V, 1);
	int avgDegree = std::count (OutDegree, OutDegree + V, avg);
	int minDegree = std::count (OutDegree, OutDegree + V, *minmax.first);
	int maxDegree = std::count (OutDegree, OutDegree + V, *minmax.second);
	std::cout << std::setprecision(1)
			  << "          Avg:  " << avg    << "\t\tOutDegree " << std::left << std::setw(14) << avg << ":  " << std::left << std::setw(14) << avgDegree << fUtil::perCent(avgDegree, V) << " %" << std::endl
			  << "     Std. Dev:  " << stdDev /*<< "\t\tOutDegree 1:  " << std::left << std::setw(14) << oneDegree  << fUtil::perCent(oneDegree, V)  << " %"*/ << std::endl
			  << "          Min:  " << *minmax.first    << "\t\tOutDegree " << std::left << std::setw(14) << *minmax.first  << ":  " << std::left << std::setw(14) << minDegree << fUtil::perCent(minDegree, V) << " %" << std::endl
			  << "          Max:  " << *minmax.second   << "\t\tOutDegree " << std::left << std::setw(14) << *minmax.second << ":  " << std::left << std::setw(14) << maxDegree << fUtil::perCent(maxDegree, V) << " %" << std::endl;
	/*if (Direction == DIRECTED)
		std::cout << "\t\t\t\t InDegree 0:  " << std::count (InDegree, InDegree + V, 0) << std::endl
				  << "\t\t\t\t InDegree 1:  " << std::count (InDegree, InDegree + V, 1) << std::endl;*/
	std::cout << std::endl;
	std::cout.imbue(std::locale());
}
