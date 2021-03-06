//#pragma once

#include "statistic.h"



/*
* Min function for percentage evaluation
*/
double min(std::vector<long double> data, int n){
	double min  = 100.0;
	for(int i = 0; i < n; i++){
		min  = (data[i] < min ? data[i] : min);
	}
	return min;
}



/*
* Max function for percentage evaluation
*/
double max(std::vector<long double> data, int n){
	double max  = 0.0;
	for(int i = 0; i < n; i++){
		max  = (data[i] > max ? data[i] : max);
	}
	return max;
}



/*
* Evaluate mean percentage of visited graph
*/
void computeMeanPercentage(std::vector<long double> Percentual, int percentCnt){
	double sum = 0;
	for(int i = 0; i < percentCnt; i++)
		sum += Percentual[i];

	printf("# Completed Visit \t: %d on %d\n", (N_TEST - percentCnt), N_TEST);
	if(percentCnt != 0)
	{
		printf("AVG Percentual \t\t: %.2f%%\n", sum / percentCnt);
		printf("MIN Percentual \t\t: %.2f%%\n", min(Percentual, percentCnt));
		printf("MAX Percentual \t\t: %.2f%%\n", max(Percentual, percentCnt));
	}
}



/*
* Evaluate Elapsed Time
*/
void computeElapsedTime(std::vector<double> topDown_times, std::vector<double> seq_times, std::vector<double> bottomUp_times){
	double sum_top = 0;
	double sum_bot = 0;
	double sum_seq = 0;

	for (int i = 0; i < N_TEST; ++i){
		sum_top += topDown_times[i];
		sum_bot += bottomUp_times[i];
		sum_seq += seq_times[i];
	}

	printf("AVG HOST TIME      \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_seq / (N_TEST), 27, 0);
	printf("AVG TOP-DOWN TIME  \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_top / (N_TEST), 27, 0);
	printf("AVG BOTTOM-UP TIME \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, sum_bot / (N_TEST), 27, 0);
	printf("AVG DEVICE TIME    \t: %c[%d;%dm%.1f%c[%dm ms\n", 27, 0, 31, (sum_bot + sum_top) / (N_TEST), 27, 0);
	printf("AVG TOTAL TIME     \t: %c[%d;%dm%.1f%c[%dm ms\n\n", 27, 0, 31, (sum_top + sum_seq + sum_bot) / (N_TEST), 27, 0);
}



/*
* Evaluate sources number at run time
*/
int EvaluateSourcesNum(float avgDeg, int N){

	int i = 0;
	float sum = 0;
	std::vector<float> FrontierStep;
	FrontierStep.push_back(1.0);

	while(FrontierStep[i] < BLOCK_FRONTIER_LIMIT)
	{
		FrontierStep.push_back( (FrontierStep[i] * 0.2) + (FrontierStep[i] * 0.8 * avgDeg) );
		i++;
	}

	for (uint k = 0; k < FrontierStep.size()-1; k++)
		sum += FrontierStep[k];

	return (int)ceil(N/sum);
}



/*
* Evaluate probability that the next node is already been visited by someone else
*/
std::vector<double> probability( int N, int Nsources, double avgDeg ){

	std::vector<double> prob;
	int PV = Nsources - 1;
	int MyPV = 1;

	prob.push_back((Nsources - 1) / (N - 1));

	while( (MyPV + PV) < N )
	{
		double Pi = (PV + (PV * avgDeg)) / (N - (MyPV * avgDeg));
		prob.push_back( Pi < 1.0 ? Pi : 1.0 );
		PV = PV + (PV * avgDeg);
		MyPV = MyPV * avgDeg;
	}

	return prob;
}