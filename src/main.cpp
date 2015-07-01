#include <iostream>
#include <datastructures.hpp>

int main() {
	dat::Matrix< complex float> M1 = eye< complex float>(5);
	dat::Matrix< complex float> M2 = (3 + I) * M1;

	std::cout << "M1: \n";
	std::cout << M1;
	std::cout << "M2: \n";
	std::cout << M2;
	std::cout << "M1+M2: \n";
	std::cout << M2 + 2 * M1 + 3 * M2;
	std::cout << "M1-M2: \n";
	std::cout << M1 - M2;
	std::cout << "M1*M2: \n";
	std::cout << M1 * M2;
	std::cout << "M2*M2: \n";
	std::cout << M2 * M2 * M2;

	std::cout << "M2[0:4,0:4] \n";
	std::cout << M2(0,4,0,4);
	std::cout << "M2[0:4,1:2] \n";
	std::cout <<  M2(0, 4, 1, 2);

	M2(0, 2, 0, 2)(0,0) = -1;
	std::cout << "M2: \n";
	std::cout<<M2;

}

