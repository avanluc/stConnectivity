//#pragma once

#include "stConn.h"


/*
* BFS on adjacency matrix performed on CPU
*/
bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target, int* Queue) {
	int left = 0, right = 1;
	std::vector<bool> Visited(nof_nodes);

	Queue[0] = source;
	Visited[source] = true;

	while (left < right)
	{
		int qNode = Queue[left++];
		for (int i = 0; i < nof_nodes; ++i)
		{
			if(adjMatrix[qNode*nof_nodes + i])
				if(i == target)
					return true;
				if (!Visited[i])
				{ 
					Visited[i] = true;
					Queue[right++] = i;
				}
		}
	}
	return (Visited[target] == true);
}



/*
* Function that choose Nsources nodes of the graph
*/
void ChooseRandomNodes(int* sources, const int V, const int Nsources, const int src, const int dst) {

	int j = 0;
	std::set<int> SourceSet;
	std::set<int>::iterator it;

	SourceSet.insert(src);

	if (Nsources > 1 )	SourceSet.insert(dst);
	if (Nsources > 2 )
	{
		while(SourceSet.size() < Nsources)
			SourceSet.insert(rand() % V);

		for (it=SourceSet.begin(); it!=SourceSet.end(); ++it, ++j)
			sources[j] = *it;
	}
	return;
}
