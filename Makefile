all:
	nvcc -arch=sm_35 -I./src/cuda  src/cuda/stConnectivity.cu --disable-warnings -o ./build/cuda/stConnectivity

seq:
	g++ -o2 src/seq/stConnectivity.cpp src/seq/stConnectivity.hpp -o ./build/seq/stConnectivity

clean:
	rm -rf ./build/*
