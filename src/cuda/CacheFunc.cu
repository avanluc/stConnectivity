#pragma once
#include <cub.cuh>



template<cub::CacheLoadModifier LDD, cub::CacheStoreModifier STD, typename T>
__device__ __forceinline__  void markWrite(T* devArray, const int index) {
	char* address = reinterpret_cast< char*>(devArray) + (index >> 3);
	const char value = cub::ThreadLoad<LDD>(address) | (1 << (7 - (index & 7)));
	cub::ThreadStore<STD>(address, value);
}



__device__ __forceinline__  void AtomicMarkWrite(int* devArray, const int index) {
	int* address = devArray + (index >> 5);
	atomicCAS(address, address[0], address[0] | (1 << (31 - (index & 31))));
}



__device__ __forceinline__  bool getBit(int* devArray, const int index) {
	return reinterpret_cast< int*>(devArray)[index >> 5] & (1 << (31 - (index & 31)));
}



template<cub::CacheLoadModifier LDD, typename T>
__device__ __forceinline__ bool markAccess(T* devArray, const int index) {
	return cub::ThreadLoad<LDD>(reinterpret_cast< char*>(devArray) + (index >> 3)) & (1 << (7 - (index & 7)));
}



template<cub::CacheLoadModifier LDD, typename T>
__device__ __forceinline__ bool markAccessInt(T* devArray, const int index) {
	return cub::ThreadLoad<LDD>(reinterpret_cast< int*>(devArray) + (index >> 5)) & (1 << (31 - (index & 31)));
}



__device__ __forceinline__ bool visitAdjiacent(int index, const int* devNode, const int* devEdge, int* BitMask)
{
	const int start = devNode[index];
	const int end 	= devNode[index + 1];
	for (int k = start; k < end; ++k)
		if(markAccessInt<cub::LOAD_CS, int>(BitMask, devEdge[k]))
			return 1;
	return 0;
}