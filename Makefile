
all:
	g++ src/stConnectivity.cpp src/stConnectivity.hpp -o ./build/stConnectivity

clean:
	rm -f ./build/*