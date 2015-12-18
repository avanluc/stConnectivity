#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <stdlib.h>
#include <set>
#include "definition.cuh"
#include "readGraph.h"
#include "Timer.cuh"

template<class T>
void PrintMatrix(T* Matrix, int size){
	for (int row = 0; row < size; ++row)
	{
		std::cout << "| ";
		for (int col = 0; col < size; ++col)
			std::cout << Matrix[ (row * size) + col ] << " ";
		std::cout << "|" << std::endl;
	}
	return;
}

bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target, int* Queue);
void ChooseRandomNodes(int* sources, const int V, const int Nsources, const int src, const int dst);
bool stConnectivity(const int* Nodes, const int* Edges, const int nof_nodes, const int source, const int target);
void PrintResults(const int test, const int source, const int target, const bool connect, const float time);
