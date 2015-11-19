#pragma once
#include <cub.cuh>


/*
* SET BIT WITH CHAR ANDRESSING
*/
template<cub::CacheLoadModifier LDD, cub::CacheStoreModifier STD, typename T>
__device__ __forceinline__  void markWrite(T* devArray, const int index) {
	char* address = reinterpret_cast< char*>(devArray) + (index >> 3);
	const char value = cub::ThreadLoad<LDD>(address) | (1 << (7 - (index & 7)));
	cub::ThreadStore<STD>(address, value);
}


/*
* SET BIT WITH INT ADDRESSING AND ATOMIC OPERATION
*/
__device__ __forceinline__  void AtomicMarkWriteInt(int* devArray, const int index) {
	int* address = devArray + (index >> 5);
	atomicCAS(address, address[0], address[0] | (1 << (31 - (index & 31))));
}


/*
* ACCESS BIT WITH CHAR ADDRESSING AND WITHOUT CACHE
*/
__device__ __forceinline__  bool getBit(int* devArray, const int index) {
	return reinterpret_cast< char*>(devArray)[index >> 3] & (1 << (7 - (index & 7)));
}


/*
* ACCESS BIT WITH CHAR ADDRESSING
*/
template<cub::CacheLoadModifier LDD, typename T>
__device__ __forceinline__ bool markAccess(T* devArray, const int index) {
	return cub::ThreadLoad<LDD>(reinterpret_cast< char*>(devArray) + (index >> 3)) & (1 << (7 - (index & 7)));
}


/*
* ACCESS BIT WITH INTEGER ADDRESSING
*/
template<cub::CacheLoadModifier LDD, typename T>
__device__ __forceinline__ bool markAccessInt(T* devArray, const int index) {
	return cub::ThreadLoad<LDD>(reinterpret_cast< int*>(devArray) + (index >> 5)) & (1 << (31 - (index & 31)));
}


/*
* VISIT INDEX NODE NEIGHBORS
*/
__device__ __forceinline__ bool visitAdjiacentBit(int index, const int* devNode, const int* devEdge, int* BitMask)
{
	const int start = devNode[index];
	const int end 	= devNode[index + 1];
	for (int k = start; k < end; ++k)
		if(markAccess<cub::LOAD_CG, int>(BitMask, devEdge[k]))
			return 1;
	return 0;
}


__device__ __forceinline__ bool visitAdjiacent(int index, const int* devNode, const int* devEdge, bool* BitMask)
{
	const int start = devNode[index];
	const int end 	= devNode[index + 1];
	for (int k = start; k < end; ++k)
		if(BitMask[devEdge[k]])
			return 1;
	return 0;
}