FLAGS=-arch=sm_35 --disable-warnings --use_fast_math 
MEM_INFO=--ptxas-options=-v

all:
	nvcc $(FLAGS) -I./src/cuda  src/cuda/stConnectivity.cu -o ./build/cuda/stConnectivity

seq:
	g++ -o2 src/seq/stConnectivity.cpp src/seq/stConnectivity.hpp -o ./build/seq/stConnectivity

clean:
	rm -rf ./build/*
