#pragma once
#include <cub.cuh>


template<cub::CacheLoadModifier LDD, typename T>
__device__ __forceinline__ bool markAccess(T* devArray, const int index) {
	return cub::ThreadLoad<LDD>(reinterpret_cast< char*>(devArray) + (index >> 3)) & (1 << (8 - (index & 7)));
}



template<cub::CacheLoadModifier LDD, cub::CacheStoreModifier STD, typename T>
__device__ __forceinline__  void markWrite(T* devArray, const int index) {
    char* address = reinterpret_cast< char*>(devArray) + (index >> 3);
    const char value = cub::ThreadLoad<LDD>(address) | (1 << (8 - (index & 7)));
    cub::ThreadStore<STD>(address, value);
}



__device__ __forceinline__  bool getBit(int* devArray, const int index) {
    int value = reinterpret_cast<int*>(devArray)[index >> 5];
    return value & (1 << (index & 31));
}



__device__ __forceinline__  void setBit(int* devArray, const int index) {
    int value = reinterpret_cast<int*>(devArray)[index >> 5];
    devArray[index >> 5] =  value | (1 << (index & 31));
}


__device__ __forceinline__ bool visitAdjiacent(int index, const int* devNode, const int* devEdge, int* BitMask, const int N)
{
	const int start = devNode[index];
	const int end = devNode[index + 1];
	for (int k = start; k < end; ++k)	
		if(markAccess<cub::LOAD_CS>(BitMask, devEdge[k]))
			return 1;
	return 0;
}