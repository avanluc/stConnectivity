#pragma once

#include <vector_types.h>	//int2

class Graph {
	public:
		int N;
		int E;
		int2* Edges;

		Graph(const int N, const int E); 
};