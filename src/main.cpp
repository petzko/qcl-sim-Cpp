#include <matrix.hpp>
#include <iostream>

using namespace std;
using namespace dat;

int main() {
	dat::Matrix<  float> M1 = eye<  float>(5);


	cout << M1;

}

