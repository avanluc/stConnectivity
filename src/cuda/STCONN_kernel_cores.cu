#pragma once

#include <assert.h>
#include "GlobalSync.cu"
#include "GlobalWrite.cu"
#include "CacheFunc.cu"
#include <curand_kernel.h>

/*
*
*/
__device__ __forceinline__ int visit(const int* devNode,
									const int*  devEdge,
									bool* BitMask,
									const int index)
{
	const int start = devNode[index];
	const int end = devNode[index + 1];
	for (int k = start; k < end; k++)
		if(BitMask[devEdge[k]] == 1)
			return 1;
	
	return 0;
}



/*
*
*/
__device__ __forceinline__ void findConnection(const int dest,
										const int index,
										const int2 current,
										const int2 destination,
										int2* devDistance,
										bool* BitMask,
										bool* Matrix,
										int& founds)
{
	if(BitMask[dest] == 1 )
	{
		if( BitMask[index] == 0 )
		{
			devDistance[index].x = destination.x + 1;
			devDistance[index].y = destination.y;
			BitMask[index] = 1;
			founds++;
		}
		else if( current.y != destination.y && destination.y < MAX_SIZE && current.y < MAX_SIZE)
		{
			Matrix[ (current.y 	   * MAX_SIZE) + destination.y ] = true;
			Matrix[ (destination.y * MAX_SIZE) + current.y     ] = true;
		}
	}
}