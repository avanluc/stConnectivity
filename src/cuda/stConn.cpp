//#pragma once

#include "stConn.h"


/*
* BFS on adjacency matrix performed on CPU
*/
bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target, int* Queue)
{
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
			{
				if(i == target)
					return true;
				if (!Visited[i])
				{ 
					Visited[i] = true;
					Queue[right++] = i;
				}
			}
		}
	}
	return (Visited[target] == true);
}



/*
* Function that choose Nsources nodes of the graph
*/
void ChooseRandomNodes(int* sources, const int V, const int Nsources, const int src, const int dst)
{

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



bool stConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target)
{
	int left = 0, right = 1;			
	int* Queue = new int[nof_nodes];
	std::vector<bool> Visited(nof_nodes);
	
	Queue[0] = source;
	Visited[source] = true;
	
	while (left < right)
	{
		int qNode = Queue[left++];

		for (int i = Nodes[qNode]; i < Nodes[qNode + 1]; ++i)
		{
			int dest = Edges[i];

			if (dest == target)
				return true;

			if (!Visited[dest]){
				Visited[dest] = true;
				Queue[right++] = dest;
			}
		}
	}
	delete Queue;
	return false;
}


void PrintResults(const int test, const int source, const int target, const bool connect, const float time)
{
	if(time < 0.0)
		printf("    \tsource: %d     \ttarget: %d      \tresult: %c[%d;%dm%s%c[%dm\n\n", 
														source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
														27, 0);
	else
		printf("#%d:\tsource: %d     \ttarget: %d      \tresult: %c[%d;%dm%s%c[%dm\t\ttime = %c[%d;%dm%.1f%c[%dm ms\n", 
														test+1, source, target, 27, 0, 31 + connect,(connect ? "true" : "false"), 
														27, 0, 27, 0, 31, time, 27, 0);
}