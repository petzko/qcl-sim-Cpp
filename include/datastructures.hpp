#ifndef _MATRIX_
#define _MATRIX_

#include <complex.h>
#include <typeinfo>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

#define DATASTRUCT_OUT_WRN(str,FILE,LINE){\
	std::cout<<"\nWARNING: (in FILE:" << FILE << "; LINE:" << LINE << ") :\n" << str << "\n";}

#define DATASTRUCT_OUT_ERR(str,FILE,LINE){\
	std::cout<<"FILE:" << FILE << " LINE:" << LINE << ": \n ERROR:" << str << "\n";}

#define explicit_cast(type,data){*((type*)&data)}
//predefine the BLAS supported datatypes -> single and double precision real and complex numbers...

const int FLT = 0;
const int DBL = 1;
const int CPLXFLT = 2;
const int CPLXDBL = 3;

std::ostream& operator<<(std::ostream& out, complex double nr);
std::ostream& operator<<(std::ostream& out, complex float nr);

template<typename _Tp>
int gettype() {

	double dbl = 0.;
	float flt = 0.;
	complex double cdbl = 0. + I;
	complex float cflt = 0. + I;

	if (typeid(_Tp) == typeid(flt))
		return FLT;
	if (typeid(_Tp) == typeid(dbl))
		return DBL;
	if (typeid(_Tp) == typeid(cflt))
		return CPLXFLT;
	if (typeid(_Tp) == typeid(cdbl))
		return CPLXDBL;

	return -1;

}

/***************************************************************
 * Usual 2D matrix of size dimI x dimJ
 ***************************************************************/
namespace dat {

template<typename _Tp>
class Matrix {

	//define the data container... -> column wise ordered
	// M*N matrix of complex numbers
	int dimI, dimJ, bytes;
	bool view;
	_Tp** data;

private:
	void init(int dim_i, int dim_j);
	/*
	 * This function returns
	 * a matrix view, or a submatrix, of the matrix m.
	 *
	 * UL element of the submatrix is (i1,j1) whereas the lower right LR is (i2,j2).
	 *
	 * The user must make sure that the following relations hold:
	 *
	 * 0<=i1<=i2<=dimI and 0<=j1<=j2M=dimJ
	 *
	 *
	 * The submatrix has n1 rows and n2 columns.
	 * where n1 = i2-i1+1; n2 = j2-j1+1;
	 *
	 * The physical number of columns in memory given by dimJ is unchanged.
	 * Mathematically, the (i,j)-th element of the new matrix is given by,
	 * m'(i,j) = m->data[(i1*dimJ + j1) + i*dimJ + j]
	 * where dimJ is the number of colums in the original matrix m!!!
	 *
	 */

	Matrix(Matrix<_Tp> m, int i1, int i2, int j1, int j2);


public:
	int getDim_i() const;
	int getDim_j() const;

	_Tp** const getMtxData() const;

	int getBytesize() const {
		return bytes;
	}

	bool isview() const {
		return view;
	}

	bool square() const {
		return dimI == dimJ;
	}

	// initialize an empty complex Nr matrix
	Matrix(int dim_i, int dim_j);
	// copy constructor
	Matrix(const Matrix<_Tp>& arg);



	~Matrix();

	/**
	 * Overload the assignment operator!! Extremely important in order to be able to correctly execute arithmetic operations on matrices with  syntax!!!
	 *
	 *
	 */
	void operator=(const Matrix<_Tp>& arg);

	/**
	 * Return a reference (implicit pointer) to the i,j-th data element
	 * so that the user can actually overwrite the corresponding entry
	 */

	_Tp& operator()(int i, int j) const;

	/**
	 * Overload the function call operator ->  indices
	 *
	 *	@param i1 - the row index of the upper left element  of the	submatrix
	 *  @param j1 - the col index of the upper left element  of the	submatrix
	 *  @param i2 - the row index of the lower right element of the submatrix;
	 *  @param j2 - the col index of the lower right element of the submatrix;
	 *  @return - a submatrix object of the current matrix with indices (i1,j1)-->(i2,j2)!
	 *
	 */
	Matrix<_Tp> operator()(int i1, int i2,int j1,int j2);

	// overload the multiplication by matrix operator...
	Matrix<_Tp> operator*(const Matrix<_Tp>& arg) const;

	//overload the addition operator
	Matrix<_Tp> operator+(const Matrix<_Tp>& arg) const;

	//overload the subtraction operator
	Matrix<_Tp> operator-(const Matrix<_Tp>& arg);

	//overload the print operator
	friend std::ostream& operator<<(std::ostream& os,
			const dat::Matrix<_Tp>& obj) {
		int dimROW = obj.getDim_i();
		int dimCOL = obj.getDim_j();
		for (int i = 0; i < dimROW; i++) {
			for (int j = 0; j < dimCOL; j++)
				os << obj.getMtxData()[i][j] << " ";
			os << "\n";
		}

		return os;
	}

	friend Matrix<_Tp> operator*(_Tp lhs, const Matrix<_Tp>& rhs) {

		dat::Matrix<_Tp> res = rhs;
		int M = rhs.getDim_i();
		int N = rhs.getDim_j();
		int type = gettype<_Tp>();
//		std::cout<< "RHS is a " <<M <<"x" << N <<" dimensional mtx. Storage space (bytes) = " << rhs.bytes <<"\n";

		//because of the templated nature
		switch (type) {
		case FLT:
			cblas_sscal(N * M, explicit_cast(float,lhs), (float*) *res.data, 1);
			break;
		case DBL:
			cblas_dscal(N * M, explicit_cast(double,lhs), (double*) *res.data,
					1);
			break;
		case CPLXFLT:
			cblas_cscal(N * M, (void*) &lhs, (void*) (*res.data), 1);
			break;
		case CPLXDBL:
			cblas_zscal(N * M, (void*) &lhs, (void*) (*res.data), 1);
			break;
		default:
			throw std::domain_error("Unsupported matrix scale operation");
		}
		return res;

	}

	friend Matrix<_Tp> operator*(Matrix<_Tp> lhs, _Tp rhs) {
		return rhs * lhs;
	}

}
;
}

template<typename _Tp>
dat::Matrix<_Tp> eye(unsigned int M) {

	dat::Matrix<_Tp> res(M, M);
	for (int i = 0; i < M; i++)
		res(i, i) = 1.0;

	return res;
}

#endif
