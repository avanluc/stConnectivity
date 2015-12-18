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


int EvaluateSourcesNum(float avgDeg, int N);
void computeElapsedTime(std::vector<double> topDown_times, std::vector<double> seq_times, std::vector<double> bottomUp_times);
void computeMeanPercentage(std::vector<long double> Percentual, int percentCnt);
double min(std::vector<long double> data, int n);
double max(std::vector<long double> data, int n);
std::vector< double > probability( int N, int Nsources, double avgDeg );
