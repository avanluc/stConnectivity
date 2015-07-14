#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <utility>
#include <climits>
#include <vector>
#include <math.h> 

using namespace std;



int main(int argc, char *argv[]){

	// Read graph from file
	int x,y;
	int  N, E;
	ifstream in (argv[1]);
	ofstream out (argv[2]);
	if(in.is_open())
	{
		// Read graph parameters and allocate memory
		in >> N >> E;
		out << N << " " <<E << endl;

		// Read graph
		for(int i = 0; i < E; i++)
		{	
			in >> x >> y;
			out << x-1 << " " << y-1 << endl;
		}
	}
	else
	{
		printf("No input file\n");
		return -1;
	}
	in.close();
	out.close();
}	