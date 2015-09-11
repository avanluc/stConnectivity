#pragma once

#include "stConn.h"


/*
* BFS on adjacency matrix performed on CPU
*/
bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target, int* Queue) {
	int left = 0, right = 1;			// left = # nodi da cui ho fatto la visita, right = # nodi incontrati durante la visita 
	std::vector<bool> Visited(nof_nodes);	// indica se un nodo è stato visitato

	Queue[0] = source;					// Si parte dalla sorgente
	Visited[source] = true;				// si marca come visitata

	while (left < right)				// fino a che ho ancora nodi dai quali non ho fatto partire la visita 
	{
		int qNode = Queue[left++];		// qNode = prossimo nodo da cui far partire la visita

		for (int i = 0; i < nof_nodes; ++i) 		// per ogni nodo nella lista di adiacenza di qNode
		{
			if(adjMatrix[qNode*nof_nodes + i])
				if (!Visited[i])					// se dest non è ancora stato visitato
				{ 
					Visited[i] = true;				// lo marco come visitato
					Queue[right++] = i;				// lo aggiungo alla coda dei nodi incontrati durante la visita 
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
