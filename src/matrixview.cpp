/*
 * matrixview.cpp
 *
 *  Created on: Jul 10, 2015
 *      Author: petzko
 */

#include <matrixview.hpp>

/**
 * #################################
 * # Matrix view Implementation!!! #
 * #################################
 */
template class dat::MatrixView<float>;
template class dat::MatrixView< complex float>;
template class dat::MatrixView<double>;
template class dat::MatrixView<complex double>;

using namespace dat;

template<class _Tp>
MatrixView<_Tp>::~MatrixView() {

}
template<class _Tp>

MatrixView<_Tp>::MatrixView(Matrix<_Tp>* parent, unsigned int id,
		unsigned int i1, unsigned int i2, unsigned int j1, unsigned int j2) {

	_dimI = i2 - i1 + 1;
	_dimJ = j2 - j1 + 1;

	_id = id;

	_viewdata = (_Tp**) malloc(sizeof(_Tp*) * _dimI);

	_parent = parent;

	if (_parent != NULL) {
		// simply store the pointers of the original matrix into view data !
		for (int i = 0; i < _dimI; i++)
			_viewdata[i] = j1 + parent->getMtxData()[i1 + i];

	} else {
		throw std::runtime_error("Cannot initialize a view from a NULL initial matrix. Please try again! \n");
	}


}
