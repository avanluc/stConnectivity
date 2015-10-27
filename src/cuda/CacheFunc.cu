#pragma once
#include <cub.cuh>

template<cub::CacheLoadModifier LDD, typename T>
__device__ __forceinline__  bool markAccess(T* devArray, const int index) {
    return cub::ThreadLoad<LDD>(reinterpret_cast<char*>(devArray) + (index >> 3)) & (1 << (index & 7));
}

/*
* (1 << (index & 7)) = 2^index a 8 bit
*/
template<cub::CacheLoadModifier LDD, cub::CacheStoreModifier STD, typename T>
__device__ __forceinline__  void markWrite(T* devArray, const int index) {
    //const char* address = reinterpret_cast< char*>(devArray) + (index >> 3);
    const char value = cub::ThreadLoad<LDD>(reinterpret_cast< char*>(devArray) + (index >> 3)) | (1 << (index & 7));
    cub::ThreadStore<STD>(reinterpret_cast< char*>(devArray) + (index >> 3), value);
}


/*
int* devArray = InputArray + Tid * Bid * 8 + Tid * 4;
int Queue[8];

reinterpret_cast<int4*>(Queue)[0] = reinterpret_cast<int4*>( devArray )[0];
reinterpret_cast<int4*>(Queue)[4] = __ldg(&reinterpret_cast<int4*>( devArray )[BLOCKDIM]);

*/