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

template<typename _Tp>
char* bytewisecopy(_Tp* val, int type, int length) {
	unsigned int bytesize = 0;

	if (type == FLT)
		bytesize = sizeof(float);
	else if (type == DBL)
		bytesize = sizeof(double);
	else if (type == CPLXFLT)
		bytesize = sizeof(complex float);
	else if (type == CPLXDBL)
		bytesize = sizeof(complex double);

	char* result = NULL;
	if (bytesize > 0) {
		char* vval = (char*) val;
		//copy only the relevant bytes!!!
		result = (char*) malloc(length * bytesize);

	}
	return result;

}

/***************************************************************
 * Usual 2D matrix of size dimI x dimJ
 *
 ***************************************************************/
namespace dat {

template<typename _Tp>
class Matrix {

	//define the data container... -> column wise ordered
	// M*N matrix of complex numbers
	int dimI, dimJ,bytes;
	_Tp* data;
public:

	int getDim_i() const {
		return dimI;
	}
	int getDim_j() const {
		return dimJ;
	}

	bool square() const {
		return dimI == dimJ;
	}

	_Tp* getMtxData() const {
		return data;
	}

	// initialize an empty complex Nr matrix
	Matrix(int dim_i, int dim_j);
	// copy constructor
	Matrix(const Matrix<_Tp>& arg);
	~Matrix();

	// overload some operators for convenience

	// make the asignment operator copy the data...
	// coppies the data from lhs to rhs...
	Matrix<_Tp> operator=(Matrix<_Tp>& arg);

	// make the function call operator retrive the i,j th data element
	_Tp& operator()(int i, int j) const;

	// overload the multiplication by matrix operator...
	Matrix<_Tp> operator*(Matrix<_Tp>& arg) const;
	//overload the addition operator
	Matrix<_Tp> operator+(Matrix<_Tp>& arg) const;

	Matrix<_Tp> operator*(_Tp lhs);
	//overload the subtraction operator
	Matrix<_Tp> operator-(Matrix<_Tp>& arg);

	//overload the print operator
	friend std::ostream& operator<<(std::ostream& os,
			const dat::Matrix<_Tp>& obj) {
		int dimROW = obj.getDim_i();
		int dimCOL = obj.getDim_j();
		for (int i = 0; i < dimROW; i++) {
			for (int j = 0; j < dimCOL; j++)
				os << obj.data[i * dimCOL + j] << " ";
			os << "\n";
		}

		return os;
	}

	friend Matrix<_Tp> operator*(_Tp lhs, Matrix<_Tp>& rhs) {

		dat::Matrix<_Tp> res = rhs;
		int M = rhs.dimI;
		int N = rhs.dimJ;
		int type = gettype<_Tp>();
//		std::cout<< "RHS is a " <<M <<"x" << N <<" dimensional mtx. Storage space (bytes) = " << rhs.bytes <<"\n";

		//because of the templated nature
		switch (type) {
		case FLT:
			cblas_sscal(N * M, explicit_cast(float,lhs), (float*) res.data, 1);
			break;
		case DBL:
			cblas_dscal(N * M, explicit_cast(double,lhs), (double*) res.data,
					1);
			break;
		case CPLXFLT:
			cblas_cscal(N * M, (void*) &lhs,  (void*)(res.data), 1);
			break;
		case CPLXDBL:
			cblas_zscal(N * M, (void*) &lhs,  (void*)(res.data), 1);
			break;
		default:
			throw std::domain_error("Unsupported matrix scale operation");
		}
		return res;

	}

	friend Matrix<_Tp> operator*(Matrix<_Tp>& lhs, _Tp rhs) {
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

/**
 template<typename _Tp>
 dat::Matrix<_Tp> scalarmultiply(dat::Matrix<_Tp> A, dat::Matrix<_Tp> B){

 int M = A.getDim_i();
 int N = A.getDim_j();
 int K = B.getDim_i();
 int L = B.getDim_j();

 if (M != K || N != L){
 throw std::length_error("Matrix Dimensions do not aggree.");
 }

 dat::Matrix<_Tp> res(M,N);

 for (int i = 0; i < M ; i++)
 for(int j = 0; j< N; j++)
 res(i,j) = A(i,j)*B(i,j);

 return res;
 }
 */

// here the dereferencing operator shall mean conjugate transpose...
/*
 template<typename _Tp>
 dat::Matrix<_Tp> operator!(dat::Matrix<_Tp>& arg) {

 if (!arg.square()) {
 std::cout
 << "Cannot transpose a non-square matrix. Aborting \n";
 throw std::domain_error("Cannot transp a non-square matrix. Abort or handle exception!");
 }

 int M = arg.getDim_i();
 dat::Matrix<_Tp> res(M, M);
 // PETZ : TODO -> use a built in GSL routine for matrix hermitian !
 for (int i = 0; i < M; i++)
 for (int j = i + 1; j < M; j++) {
 res(j, i) = (_Tp) conj(arg(i, j));
 res(i, j) = (_Tp) conj(arg(j, i));
 }
 return res;

 }
 */
/****************************************************************
 * A 3D matrix of size dimI x dimJ x dimK (usually representing a 3D grid)
 *
 ***************************************************************/

/***
 namespace dat {

 template<typename _Tp>
 class Grid {

 //define the data container... -> column wise ordered
 // M*N matrix of complex numbers
 std::vector<int> dim;
 std::vector<int> offset;
 _Tp* data;

 public:

 inline int getDim(int i) const {
 return dim[i];
 }
 _Tp* getMtxData() const {
 return data;
 }

 _Tp* getVector(int direction, int pos_1, int pos_2) const;

 // initialize an empty complex Nr matrix
 Grid(int dim_i, int dim_j, int dim_k);
 Grid(int dim_i, int dim_j, int dim_k,_Tp val);
 // copy constructor
 Grid(const Grid<_Tp>& arg);
 ~Grid();

 // overload some operators for convenience

 // make the asignment operator copy the data...
 // coppies the data from lhs to rhs...
 Grid<_Tp>& operator=(Grid<_Tp> arg);

 // make the function call operator retrive the i,j th data element
 _Tp& operator()(int i, int j, int k) const;

 //overload the addition operator
 Grid<_Tp> operator+(Grid<_Tp> arg) const;

 //overload the subtraction operator
 Grid<_Tp> operator-(Grid<_Tp> arg) const;

 Grid<_Tp> operator*(Grid<_Tp> arg) const;

 // overload the multiplication operator
 //	Grid<_Tp> operator*(Grid<_Tp> arg) const;

 //overload the print operator

 friend std::ostream& operator<<(std::ostream& os,
 const dat::Grid<_Tp>& obj) {

 int dimI = obj.getDim(0);
 int dimJ = obj.getDim(1);
 int dimK = obj.getDim(2);

 for (int i = 0; i < dimI; i++) {
 for (int j = 0; j < dimJ; j++) {
 for (int k = 0; k < dimK; k++) {
 os << obj.data[k + j * dimI + i * dimJ * dimI] << " ";
 }

 os << "\n";
 }

 os << "\n";
 }

 return os;
 }

 }
 ;

 }

 template<typename _Tp>
 dat::Grid<_Tp> operator*(const _Tp & lhs, dat::Grid<_Tp>& rhs) {

 dat::Grid<_Tp> res = rhs;
 int dimI = res.getDim(0);
 int dimJ = res.getDim(1);
 int dimK = res.getDim(2);

 for (int i = 0; i < dimI; i++)
 for (int j = 0; j < dimJ; j++)
 for(int k = 0;k<dimK;k++)
 res(i, j,k) *= lhs;

 return res;

 }

 template<typename _Tp>
 dat::Grid<_Tp> operator*(dat::Grid<_Tp>& lhs, const _Tp& rhs) {
 return rhs * lhs;
 }
 */
#endif
