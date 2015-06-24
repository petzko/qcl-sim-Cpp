#include <iostream>
#include <datastructures.hpp>

std::ostream& operator<<(std::ostream& out, complex double nr) {

	out << "(" << creal(nr) << ";" << cimag(nr) << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, complex float nr) {
	out << "(" << creal(nr) << ";" << cimag(nr) << ")";
	return out;
}

int main() {

	std::cout << "hello world! \n";

	dat::Matrix< complex double> M1 = eye< complex double>(2);
	dat::Matrix< complex double> M2 = (3 +I)* M1;

	std::cout<< "M1: \n"; std::cout<<M1;
	std::cout<< "M2: \n";	std::cout << M2;
	std::cout<< "M1+M2: \n"; std::cout<< M1+M2 ;
	std::cout<< "M1-M2: \n"; std::cout<< M1-M2 ;
	std::cout<< "M1*M2: \n"; std::cout<< M1*M2;
	std::cout<< "M2*M2: \n"; std::cout<< M2*M2;




}

