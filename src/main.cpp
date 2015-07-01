#include <iostream>
#include <datastructures.hpp>

using namespace std;
using namespace dat;

int main() {
	dat::Matrix<  float> M1 = eye<  float>(5);

	dat::Matrix< float> m1view = M1(0,1,0,1);
	m1view(0,1) = 1; //<- shall change the state of M1

	cout << M1;

	cout << m1view;

	cout << M1;

}

