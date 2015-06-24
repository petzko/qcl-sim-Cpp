#include <iostream>
#include <datastructures.hpp>


int main() {

	std::cout << "hello world! \n";

	dat::Matrix< complex float> M1 = eye< complex float>(10);
	dat::Matrix< complex float> M2 = (3 +I)* M1;

	std::cout<< "M1: \n"; std::cout<<M1;
	std::cout<< "M2: \n";	std::cout << M2;
	std::cout<< "M1+M2: \n"; std::cout<< M1+M2 ;
	std::cout<< "M1-M2: \n"; std::cout<< M1-M2 ;
	std::cout<< "M1*M2: \n"; std::cout<< M1*M2;
	std::cout<< "M2*M2: \n"; std::cout<< M2*M2;




}

