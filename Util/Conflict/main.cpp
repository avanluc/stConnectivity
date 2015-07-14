#include <iostream>
#include <algorithm>

int main() {
	int Array[32];
	for (int i = 1; i <= 32; i++) {
		
		std::fill(Array, Array + 32, 0);
		
		for (int j = 0; j < 32; j++)
			Array[ (i * j) % 32 ]++; 
	
		int max = -1;
		for (int j = 0; j < 32; j++)
			max = std::max(max, Array[ j ]);
		std::cout << i << "  max Conflict: " << max << std::endl;
	}	
}
