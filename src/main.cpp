#include <iostream>
#include <datastructures.hpp>



int main(){

	std::cout<<"hello world! \n" ; 

	dat::Matrix< double> M1 = eye<double>(5);
	dat::Matrix< double> M2 = 3.*M1;

	std::cout<< "M1: \n"; std::cout<<M1; 
	std::cout<< "M2: \n";	std::cout << M2; 
	std::cout<< "M1+M2: \n"; std::cout<< M1+M2 ; 
	std::cout<< "M1-M2: \n"; std::cout<< M1-M2 ; 
	std::cout<< "M1*M2: \n"; std::cout<< M1*M2;



}



