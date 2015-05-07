all:
	nvcc -arch=sm_30 -I./src/cuda  src/cuda/stConnectivity.cu -O3 --disable-warnings -o ./build/cuda/stConnectivity

seq:
	g++ -o2 src/seq/stConnectivity.cpp src/seq/stConnectivity.hpp -o ./build/seq/stConnectivity

clean:
	rm -rf ./build/*
