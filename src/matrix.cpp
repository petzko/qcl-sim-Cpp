/*
 * matrix.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: petz
 */

#include <datastructures.hpp>

using namespace dat;

/**
 * the standard constructor - takes a row and column dimension.
 */
template<typename _Tp>
Matrix<_Tp>::Matrix(int dim_i, int dim_j) {

	dimI = dim_i;
	dimJ = dim_j;
	data = (_Tp*) calloc(dimI * dimJ, sizeof(_Tp));
	bytes = dimI*dimJ*sizeof(_Tp);

}

/**
 *	Copy constructor
 */
template<typename _Tp>
Matrix<_Tp>::Matrix(const Matrix<_Tp>& arg) :
		Matrix(arg.getDim_i(), arg.getDim_j()) {


	switch (gettype<_Tp>()) {
	case FLT:
		cblas_scopy(getDim_i() * getDim_j(), (float*) arg.getMtxData(), 1,
				(float*) this->data, 1);
		break;

	case DBL:
		cblas_dcopy(dimI * dimJ, (double*) arg.getMtxData(), 1,
				(double*) this->data, 1);
		break;
	case CPLXFLT:
		cblas_ccopy(dimI * dimJ, (void*) arg.getMtxData(), 1,
				(void*) this->data, 1);

		break;
	case CPLXDBL:
		cblas_zcopy(getDim_i() * getDim_j(), (void*) arg.getMtxData(), 1,
				(void*) this->data, 1);

		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");

	}


}
/*
 * Destructor
 */
template<typename _Tp>
Matrix<_Tp>::~Matrix() {
	if (data != NULL)
		free(data);
}

// make the function call operator retrive the i,j th data element
template<typename _Tp>
_Tp& Matrix<_Tp>::operator()(int i, int j) const {

	if (i >= dimI || i < 0) {
		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot retrieve Matrix element. ",
				__FILE__, __LINE__)
		throw std::out_of_range("row index out of bounds.");

	}

	if (j >= dimJ || j < 0) {
		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot retrieve Matrix element. ",
				__FILE__, __LINE__)
		throw std::out_of_range("column index out of bounds.");
	}

	return data[i * dimJ + j];

}
/* 
 overload the assignment operator!
 *
 */
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator=(Matrix<_Tp>& arg) {

	if (dimI != arg.getDim_i() || dimJ != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	switch (gettype<_Tp>()) {
	case FLT:
		cblas_scopy(getDim_i() * getDim_j(), (float*) arg.getMtxData(), 1,
				(float*) data, 1);
		break;
	case DBL:
		cblas_dcopy(getDim_i() * getDim_j(), (double*) arg.getMtxData(), 1,
				(double*) data, 1);
		break;
	case CPLXFLT:
		cblas_ccopy(getDim_i() * getDim_j(), (void*) arg.getMtxData(), 1,
				(void*) data, 1);
		break;
	case CPLXDBL:
		cblas_zcopy(getDim_i() * getDim_j(), (void*) arg.getMtxData(), 1,
				(void*) data, 1);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");

	}

	return (*this);

}

// overload the multiplication by matrix operator...
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator*(Matrix<_Tp>& arg) const {

	if (dimJ != arg.getDim_i()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	int L = arg.getDim_j();
	Matrix<_Tp> res(dimI, L);

	complex float alpha_f = 1.0;
	complex double alpha_d = 1.0;
	complex float beta_f = 0.;
	complex double beta_d = 0.;
	switch (gettype<_Tp>()) {
	case FLT:
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimI,
				arg.getDim_j(), dimJ, 1.f, (float*) data, dimI,
				(float*) arg.getMtxData(), arg.getDim_i(), 0.f,
				(float*) res.getMtxData(), dimI);
		break;
	case DBL:
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimI,
				arg.getDim_j(), dimJ, 1., (double*) data, dimI,
				(double*) arg.getMtxData(), arg.getDim_i(), 0.,
				(double*) res.getMtxData(), dimI);
		break;
	case CPLXFLT:
		cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimI,
				arg.getDim_j(), dimJ, (void*) &alpha_f, (void*) data, dimI,
				(void*) arg.getMtxData(), arg.getDim_i(), (void*) &beta_f,
				(void*) res.getMtxData(), dimI);
		break;

	case CPLXDBL:
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimI,
				arg.getDim_j(), dimJ, (void*) &alpha_d, (void*) data, dimI,
				(void*) arg.getMtxData(), arg.getDim_i(), (void*) &beta_d,
				(void*) res.getMtxData(), dimI);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix-matrix mutiplication!");

	}

	return res;
}

//overload the addition operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator+(Matrix<_Tp>& arg) const {

	if (dimI != arg.getDim_i() || dimJ != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res = arg;

	complex double alpha_d = 1.;
	complex float alpha_f = 1.;
	int len = dimI * dimJ;

	switch (gettype<_Tp>()) {
	case FLT:
		cblas_saxpy(len, 1.0, (float*) data, 1, (float*) res.data, 1);
		break;
	case DBL:
		cblas_daxpy(len, 1.0, (double*) data, 1, (double*) res.data, 1);
		break;
	case CPLXFLT:
		cblas_caxpy(len, (void*) &alpha_f, (void*) data, 1,
				(void*) res.getMtxData(), 1);
		break;
	case CPLXDBL:
		cblas_zaxpy(len, (void*) &alpha_d, (void*) data, 1,
				(void*) res.getMtxData(), 1);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix-matrix addition!");

	}

	return res;
}

//overload the subtraction operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator-(Matrix<_Tp>& arg) {

	if (dimI != arg.getDim_i() || dimJ != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res = ((_Tp) -1.) * arg;
	return (*this) + res;

}

template class dat::Matrix<float>;
template class dat::Matrix< complex float>;
template class dat::Matrix<double>;
template class dat::Matrix<complex double>;

