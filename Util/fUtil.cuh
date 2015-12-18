#pragma once

#include <sys/syscall.h>
#include <sys/resource.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <limits>
#include <numeric>
#include <cmath>
#include <cuda_runtime.h>	// error
#include <iomanip>      	// setprecision
#include <cstdlib>			// exit

#include <locale>			// numpunct<char>
#include <ostream>			// color


#define DIV(a,b) (((a) + (b) - 1)/(b))		// q = ceil(x / y)   =>   q = (x + y - 1) / y;


// ----------- META-PROGRAMMING -------------------

template <int _N, int CURRENT_VAL = _N, int COUNTL = 0>
struct _Log2 {
	enum { VALUE = _Log2< _N, (CURRENT_VAL >> 1), COUNTL + 1>::VALUE };
};

template <int _N, int COUNTL>
struct _Log2<_N, 0, COUNTL> {
	enum { VALUE = (1 << (COUNTL - 1) < _N) ? COUNTL : COUNTL - 1 };
};


template <int POW, int COUNTL = POW>
struct _Pow2 {
	enum { VALUE = POW == 0 ? 1 : _Pow2< (POW << 1), COUNTL - 1 >::VALUE };
};

template <int POW>
struct _Pow2<POW, 0> {
	enum { VALUE = POW };
};

template <int MOD>
struct _Mod2 {
	enum { VALUE = MOD - 1 };
};

// NEAREST POW 2

template <int NUM, int C>
struct _OR_SHIFT {
	enum { VALUE = _OR_SHIFT<NUM | (NUM >> C), C * 2>::VALUE};
};

template <int NUM>
struct _OR_SHIFT<NUM, 32> {
	enum { VALUE = NUM + 1 };
};

template <int NUM>
struct _NearestPow2 {
	enum { VALUE = _OR_SHIFT<NUM - 1, 1>::VALUE };
};
	

namespace fUtil {


	struct myseps : std::numpunct<char> {
		// use space as separator
		char do_thousands_sep() const { return ','; }

		// digits are grouped by 3 digits each
		std::string do_grouping() const { return "\3"; }
	};

	void memInfoPrint(size_t total, size_t free, size_t Req);
	void cudaStatics();
	void memInfoCPU(int Req);
	void memInfoCUDA(int Req);

	std::string extractFileName(std::string s);
	std::string NumToString(int x);

	// --------------------------- MATH ---------------------------------------------------

	unsigned nearestPower2(unsigned v);
	unsigned log_2(unsigned n);
	bool isPositiveInteger(const std::string& s);
	float perCent(int part, int max);

	template<class T>
	void MinMax_Range(T* Input, const int size, T& min, T& max);
	
	template<class T>
	float average(T* Input, const int size);

	template<class T>
	float stdDeviation(T* Input, const int size);

	template<class T>
	void FisherYatesRandomize(T* array, const int N);

	void genRandomPermutation(int*& Randomize, const int N);

	// --------------------------- ARRAY ---------------------------------------------------
	template<class T>
	void printArray(T Array[], const int size, std::string text = "", bool debug = true, const char sep = ' ', T inf = std::numeric_limits<T>::max());
	
	void printArray2(int2 Array[], const int size, std::string text, bool debug = true);
	
	template<class H, class D>
	void Compare(H* HostArray, D* DevArray, const int N, bool PRINT = false);

	template<class T>
	void printMatrix(T** Matrix, const int ROW, const int COL, std::string text, bool debug, T inf);

#if defined(__NVCC__)
	inline void __getLastCudaError(const char *errorMessage, const char *file, const int line);

	template<class H, class D>
	void cudaCompare(H* hostArray, D* devArray, const int N, const int CONTROL);

	template<class T>
	void printCudaArray(T* devArray, const int size, std::string text = "", bool debug = true, const char sep = ' ', T inf = std::numeric_limits<T>::max());

	template<class T>
	__global__ void scatterKernel(const int* toScatter, const int nof_items, T* Dest, const T value);

	template<class T>
	__global__ void fillKernel(T* toFill, const int nof_items, const T value);
#endif

	// --------------------------- IMPLEMENTATION ---------------------------------------------------
	template<class T>
	void MinMax_Range(T* Input, const int size, T& min, T& max) {
		min = Input[0];
		max = Input[0];
		for (int i = 1; i < size; i++) {
			if (Input[i] < min)
				min = Input[i];
			if (Input[i] > max)
				max = Input[i];
		}
	}
	
	template<class T>
	float average(T* Input, const int size) {
		const T sum = std::accumulate(Input, Input + size, 0);
		return (float) sum / size;
	}

	template<class T>
	float stdDeviation(T* Input, const int size) {
		const float avg = average(Input, size);
		return stdDeviation(Input, size, avg);
	}
	
	template<class T>
	float stdDeviation(T* Input, const int size, float avg) {
		float sum = 0;
		for (int i = 0; i < size; ++i)
			sum += std::pow(Input[i] - avg, 2);
		return std::sqrt(sum / size);
	}

	template<class T>
	void FisherYatesRandomize(T* array, const int N) {
		for (int i = 0; i < N - 1; i++) {
			int j = i + mt::randf_co() * (N - i);
			std::swap<T>(array[i], array[j]);
		}
	}

	template<class H, class D>
	void Compare(H* HostArray, D* DevArray, const int N, bool PRINT) {
		printArray<H>(HostArray, N,   "Host Array:\t", PRINT);
		printArray<D>(DevArray, N,  "Device Array:\t", PRINT);

		for (int i = 0; i < N; i++) {
			if (HostArray[i] != (H) DevArray[i] && !(HostArray[i] == std::numeric_limits<H>::max() && DevArray[i] == std::numeric_limits<D>::max()))
				error("Different Value: index " << i << "\tHost: " << HostArray[i] << "\tDevice: " << DevArray[i]);
		}
		std::cout << "\n<> CORRECT\n";
	}

	template<class T>
	void printArray(T Array[], const int size, std::string text, bool debug, const char sep, T inf) {
		if (!debug)
			return;

		std::cout << text;
		for (int i = 0; i < size; i++)
			std::cout << ((Array[i] == inf) ? "inf" : std::string(NumToString(Array[i]))) << sep;
		std::cout << std::endl << std::endl;
	}
	

	template<class T>
	void printMatrix(T** Matrix, const int ROW, const int COL, std::string text, bool debug, T inf) {
		if (!debug)
			return;

		std::cout << text;
		for (int i = 0; i < ROW; i++)
			printArray(Matrix[i * COL], COL, "\n", true, inf, '\t');
		std::cout << std::endl << std::endl;
	}

#if defined(__NVCC__)
	inline void __getLastCudaError(const char *errorMessage, const char *file, const int line) {
		cudaError_t err = cudaGetLastError();
		if (cudaSuccess != err) {
			std::cerr << std::endl << file << "(" << line << ") : getLastCudaError() CUDA error : " << errorMessage << " : (" << (int) err << ") " << cudaGetErrorString(err) << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	template<class H, class D>
	void cudaCompare(H* hostArray, D* devArray, const int N, const int CONTROL) {
		if (!CONTROL)
			return;
		D* ArrayCMP = new D[N];
		cudaMemcpy(ArrayCMP, devArray, N * sizeof (D), cudaMemcpyDeviceToHost);
		cudaError("Copy To Host");

		Compare<H, D>(hostArray, ArrayCMP, N);
		delete[] ArrayCMP;
	}

	template<class T>
	void printCudaArray(T* devArray, const int size, std::string text, bool debug, const char sep, T inf) {
		if (!debug)
			return;

		T* hostArray = new T[size];
		cudaMemcpy(hostArray, devArray, size * sizeof (T), cudaMemcpyDeviceToHost);
		cudaError("Copy To Host");

		printArray<T>(hostArray, size, text, debug, sep, inf);
		delete hostArray;
	}

	template<class T>
	__global__ void scatterKernel(const int* toScatter, const int nof_items, T* Dest, const T value) {
		const int ID = blockIdx.x * blockDim.x + threadIdx.x;

		if (ID < nof_items)
			Dest[ toScatter[ID] ] = value;
	}

	template<class T>
	__global__ void fillKernel(T* toFill, const int nof_items, const T value) {
		const int ID = blockIdx.x * blockDim.x + threadIdx.x;

		if (ID < nof_items)
			toFill[ ID ] = value;
	}
#endif
}
