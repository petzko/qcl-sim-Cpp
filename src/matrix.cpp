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

std::ostream& operator<<(std::ostream& out, complex double nr) {

	out << "(" << creal(nr) << ";" << cimag(nr) << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, complex float nr) {
	out << "(" << creal(nr) << ";" << cimag(nr) << ")";
	return out;
}

template<typename _Tp>
int Matrix<_Tp>::getDim_i() const {
	return dimI;
}

template<typename _Tp>
int Matrix<_Tp>::getDim_j() const {
	return dimJ;
}

template<typename _Tp>
_Tp** const Matrix<_Tp>::getMtxData() const {
	return data;
}

template<typename _Tp>
void Matrix<_Tp>::init(int dim_i, int dim_j) {

	dimI = dim_i;
	dimJ = dim_j;
	view = false;

	_Tp* tmp_data = (_Tp*) calloc(dimI * dimJ, sizeof(_Tp));
	data = (_Tp**) calloc(dimI, sizeof(_Tp*));

	for (int i = 0; i < dimI; i++)
		data[i] = tmp_data + i * dimJ;

	bytes = dimI * dimJ * sizeof(_Tp);
}

template<typename _Tp>
Matrix<_Tp>::Matrix(int dim_i, int dim_j) {

	init(dim_i, dim_j);

}

template<typename _Tp>
dat::Matrix<_Tp>::Matrix(Matrix<_Tp> m, int i1, int i2, int j1, int j2) {

	view = true;
	int k1 = i1;
	int k2 = j1;
	int n1 = i2 - i1 + 1;
	int n2 = j2 - j1 + 1;

	dimI = n1;
	dimJ = n2;

	data = (_Tp**) calloc(n1, sizeof(_Tp));
	bytes = n1 * n2 * sizeof(_Tp);

	// copy the pointers NOT the data !!!!
	for (int i = 0; i < n1; i++)
		data[i] = m.getMtxData()[i + k1] + k2;

}

template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator()(int i1, int i2, int j1, int j2) {

	// index out of bounds error check
	if (i1 < 0 || i2 - i1 < 0 || i2 >= dimI || j1 < 0 || j2 - j1 < 0
			|| j2 >= dimJ) {

		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot generate submatrix-view with that dimensions. ",
				__FILE__, __LINE__)
		throw std::out_of_range("array index out of bounds.");

	}
	// calls the COPY CONSTRUCTOR!!!
	Matrix<_Tp> res = Matrix(*this, i1, i2, j1, j2);
	return res;

}

/**
 *	Copy constructor
 */
template<typename _Tp>
Matrix<_Tp>::Matrix(const Matrix<_Tp>& arg) {

	if (arg.isview()) {
		view = true;

		dimI = arg.getDim_i();
		dimJ = arg.getDim_j();

		data = (_Tp**) calloc(dimI, sizeof(_Tp));
		bytes = dimI * dimJ * sizeof(_Tp);

		for (int i = 0; i < dimI; i++)
			data[i] = arg.getMtxData()[i];
		return;
	}


	init(arg.dimI, arg.dimJ);

	switch (gettype<_Tp>()) {
	case FLT:
		for (int i = 0; i < getDim_i(); i++) {
			cblas_scopy(getDim_j(), (float*) arg.getMtxData()[i], 1,
					(float*) this->getMtxData()[i], 1);
		}
		break;

	case DBL:
		for (int i = 0; i < getDim_i(); i++) {
			cblas_dcopy(getDim_j(), (double*) arg.getMtxData()[i], 1,
					(double*) this->getMtxData()[i], 1);
		}
		break;
	case CPLXFLT:

		for (int i = 0; i < getDim_i(); i++) {
			cblas_ccopy(getDim_j(), (void*) arg.getMtxData()[i], 1,
					(void*) this->getMtxData()[i], 1);
		}
		break;
	case CPLXDBL:
		for (int i = 0; i < getDim_i(); i++) {
			cblas_zcopy(getDim_j(), (void*) arg.getMtxData()[i], 1,
					(void*) this->getMtxData()[i], 1);
		}
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

	if (view)
		// could cause a memeory leak !!! as pointers to data are not freed!
		data =  NULL;
	else {
		free(data); // free the whole data!
// no need for a for loop as data is contiguously allocated
//		for (int i = 0; i < dimI; i++)
//			free(data[i]);
//		free(data);
	}
}

// make the function call operator retrive the i,j th data element
template<typename _Tp>
_Tp& Matrix<_Tp>::operator()(int i, int j) const {

	if (i >= this->getDim_i() || i < 0) {
		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot retrieve Matrix element. ",
				__FILE__, __LINE__)
		throw std::out_of_range("row index out of bounds.");

	}

	if (j >= this->getDim_j() || j < 0) {
		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot retrieve Matrix element. ",
				__FILE__, __LINE__)
		throw std::out_of_range("column index out of bounds.");
	}

		return data[i][j];

}
/* 
 * overload the assignment operator!
 *
 */
template<typename _Tp>
void Matrix<_Tp>::operator=(const Matrix<_Tp>& arg) {

	if (arg.isview())
		for (int i = 0; i < dimI; i++){
			data[i] = arg.getMtxData()[i];
		}

	if (this->getDim_i() != arg.getDim_i()
			|| this->getDim_j() != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	switch (gettype<_Tp>()) {
	case FLT:
		for (int i = 0; i < getDim_i(); i++) {
			cblas_scopy(getDim_j(), (float*) arg.getMtxData()[i], 1,
					(float*) this->getMtxData()[i], 1);
		}
		break;

	case DBL:
		for (int i = 0; i < getDim_i(); i++) {
			cblas_dcopy(getDim_j(), (double*) arg.getMtxData()[i], 1,
					(double*) this->getMtxData()[i], 1);
		}
		break;
	case CPLXFLT:

		for (int i = 0; i < getDim_i(); i++) {
			cblas_ccopy(getDim_j(), (void*) arg.getMtxData()[i], 1,
					(void*) this->getMtxData()[i], 1);
		}
		break;
	case CPLXDBL:
		for (int i = 0; i < getDim_i(); i++) {
			cblas_zcopy(getDim_j(), (void*) arg.getMtxData()[i], 1,
					(void*) this->getMtxData()[i], 1);
		}
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");

	}
}

// overload the multiplication by matrix operator...
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator*(const Matrix<_Tp>& arg) const {

	if (this->getDim_j() != arg.getDim_i()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	int L = arg.getDim_j();
	Matrix<_Tp> res(this->getDim_i(), L);

	complex float alpha_f = 1.0;
	complex double alpha_d = 1.0;
	complex float beta_f = 0.;
	complex double beta_d = 0.;
	switch (gettype<_Tp>()) {
	case FLT:
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), 1.f,
				(float*) *this->getMtxData(), this->getDim_i(),
				(float*) *arg.getMtxData(), arg.getDim_i(), 0.f,
				(float*) *res.getMtxData(), this->getDim_i());
		break;
	case DBL:
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), 1.,
				(double*) *this->getMtxData(), this->getDim_i(),
				(double*) *arg.getMtxData(), arg.getDim_i(), 0.,
				(double*) *res.getMtxData(), this->getDim_i());
		break;
	case CPLXFLT:
		cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), (void*) &alpha_f,
				(void*) *this->getMtxData(), this->getDim_i(),
				(void*) *arg.getMtxData(), arg.getDim_i(), (void*) &beta_f,
				(void*) *res.getMtxData(), this->getDim_i());
		break;

	case CPLXDBL:
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), (void*) &alpha_d,
				(void*) *this->getMtxData(), this->getDim_i(),
				(void*) *arg.getMtxData(), arg.getDim_i(), (void*) &beta_d,
				(void*) *res.getMtxData(), this->getDim_i());
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix-matrix mutiplication!");

	}

	return res;
}

//overload the addition operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator+(const Matrix<_Tp>& arg) const {

	if (this->getDim_i() != arg.getDim_i()
			|| this->getDim_j() != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res = arg;

	complex double alpha_d = 1.;
	complex float alpha_f = 1.;

	switch (gettype<_Tp>()) {
	case FLT:
		for (int i = 0; i < getDim_i(); i++)
			cblas_saxpy(this->getDim_j(), 1.0, (float*) this->getMtxData()[i],
					1, (float*) res.getMtxData()[i], 1);

		break;
	case DBL:
		for (int i = 0; i < getDim_i(); i++)
			cblas_daxpy(this->getDim_j(), 1.0, (double*) this->getMtxData()[i],
					1, (double*) res.getMtxData()[i], 1);
		break;
	case CPLXFLT:
		for (int i = 0; i < getDim_i(); i++)
			cblas_caxpy(this->getDim_j(), (void*) &alpha_f,
					(void*) this->getMtxData()[i], 1,
					(void*) res.getMtxData()[i], 1);
		break;
	case CPLXDBL:
		for (int i = 0; i < getDim_i(); i++)
			cblas_zaxpy(this->getDim_j(), (void*) &alpha_d,
					(void*) this->getMtxData()[i], 1,
					(void*) res.getMtxData()[i], 1);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix-matrix addition!");

	}

	return res;

}

//overload the subtraction operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator-(const Matrix<_Tp>& arg) {

	if (this->getDim_i() != arg.getDim_i()
			|| this->getDim_j() != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res = ((_Tp) -1.) * arg;
	return (*this) + res;

}

template class dat::Matrix<float>;
template class dat::Matrix< complex float>;
template class dat::Matrix<double>;
template class dat::Matrix<complex double>;

