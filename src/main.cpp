#include <iostream>
#include <datastructures.hpp>

void testview1(dat::Matrix<complex float> M){
	dat::Matrix<complex float> view = M(0,0,1,1);
	view(0,0) = 10;
	std::cout<<view;
}

dat::Matrix<complex float>* testview2(){

	dat::Matrix<complex float> Mtx = eye< complex float>(5);
	dat::Matrix<complex float>* view = new dat::Matrix<complex float >(Mtx,0,0,3,3);
	return view;
}


int main() {
	dat::Matrix< complex float> M1 = eye< complex float>(5);
	dat::Matrix< float> * M2 = new dat::Matrix<float>(1,2);



	M1(0,0) = 0;
	M1(1,1) = 1;
	M1(2,2) = 2;
	M1(3,3) = 3;
	M1(4,4) = 4;


//	dat::Matrix< complex float> M2 = 3*M1; //   (3.+3.*I)*M1;
//	std::cout<< "M1: \n"; std::cout<<M1;
//	std::cout<< "M2: \n";	std::cout << M2;
//	std::cout<< "M1+M2: \n"; std::cout<< M2+2*M1 +3*M2 ;
//	std::cout<< "M1-M2: \n"; std::cout<< M1-M2 ;
//	std::cout<< "M1*M2: \n"; std::cout<< M1*M2;
//	std::cout<< "M2*M2: \n"; std::cout<< M2*M2*M2;
//

	std::cout<<"before test view: \n";
	std::cout<< M1;
	std::cout<<"*******\n";
	testview1(M1);

	std::cout<<"*******\n";
	std::cout<<"after test view: \n";
	std::cout<< M1;

	dat::Matrix<complex float>* view2= testview2();
	std::cout<<"view2: \n";
	std::cout<<*view2;

}

