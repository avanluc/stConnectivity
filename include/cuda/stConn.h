#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <ctime>
#include "definition.cuh"
#include "readGraph.h"
#include "Timer.cuh"

// Edge structure
struct edge{
	int x;
	int y;
};

// Stuct for ordered nodes
struct ONodes{
	int id;
	int degree;
};


int compare(const void *x, const void *y);
int ONodesCompare(const ONodes &a, const ONodes &b);
int EvaluateSourcesNum(float avgDeg, int N);
bool MatrixBFS(const bool* adjMatrix, const int nof_nodes, const int source, const int target, int* Queue);
void ChooseNodes(int* sources, std::vector< ONodes > OrderedNodes, const int nof_distNodes, const int source, const int target);
void ChooseRandomNodes(int* sources, const int* Nodes, const int nof_nodes, const int nof_distNodes, const int source, const int target);
void computeElapsedTime(std::vector<double> par_times, std::vector<double> seq_times, int connectCnt);
void computeMeanPercentage(std::vector<long double> Percentual, int percentCnt);
void GraphToCSR(const edge* graph, int* vertex, int* edges, int N, int E);
void ReadGraph(char *file, edge *graph, int *N, int *E);
double min(std::vector<long double> data, int n);
double max(std::vector<long double> data, int n);
std::vector< ONodes > OrderNodes(const int* Nodes, const int nof_nodes);