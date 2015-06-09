#include "graph.hpp"

Graph::Graph(const int N, const int E){
	this->N = N;
	this->E = E;
	Edges = new int2[E];	
}