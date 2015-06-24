#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>
#include <stdexcept>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <stdlib.h>
#include <complex.h>



#define MATRIX_OUT_WRN(str,FILE,LINE){\
	std::cout<<"\nWARNING: (in FILE:" << FILE << "; LINE:" << LINE << ") :\n" << str << "\n";}
#define MATRIX_OUT_ERR(str,FILE,LINE){\
	std::cout<<"\n ERROR: (in FILE:" << FILE << "; LINE:" << LINE << ") :\n" << str << "\n";}

/**
* override the << operator for complex double variables! print the complex number z = a+ib  as a touple (a,b)
* this function MUST be outside the namespace mtx!
*/

std::ostream& operator<<(std::ostream& ostr,  complex double z){
	ostr <<"("<<creal(z) << ";" <<cimag(z) << ")"; return ostr; }

namespace mtx{
	/**
	 * cplxMatrix<M,N> defines a M(rows) x N(cols), where M,N >= 1,  matrix of complex numbers.
	 * The class overrieds basic matrix operation as matrix matrix and matrix vector add,substract,multiply, scale, scalar product etc, as well as
	 * standard routines for direct acces to matrix elementsi (real and imaginary part), printing and others.
	 * The implementation serves as a wrapper to the more low-level gsl (GNU Scientific Library) gsl_vector_complex, element.
	 *
	 **/

	class cplxMatrix{


		//define the data container... -> row major M*N matrix of complex numbers
		complex double * data = NULL;
		size_t M ; size_t N;
		public:
		// initialize an empty complex Nr matrix and set to zero
		cplxMatrix(size_t M , size_t N)
		{
			if ( M*N <= 0 )
			{	MATRIX_OUT_ERR("Incorrect matrix dimension. Aborting",__FILE__,__LINE__)
				throw std::invalid_argument("Matrix initialization failed!Incorrenct matrix dimensions!");
			}
			this->M = M ; this->N = N ;
			this->data = (complex double*) malloc(M*N*sizeof(complex double));

		}

		cplxMatrix(const cplxMatrix& obj): cplxMatrix(obj.rows(),obj.cols()) {


			  for(size_t i = 0 ; i < this->M ; i++)
	                        for (size_t j = 0 ; j < this->N  ;j++)
        	                        data[i*N+j] = obj.data[i*N+j];

		}


		// destructor -> call the gsl vector complex free impl.
		~cplxMatrix(){ if(data != NULL) { free(data) ; data = NULL;} }


		// overload some operators for convenience
		// make the asignment operator copy the data...
		// coppies the data from lhs to rhs...

		cplxMatrix& operator= (cplxMatrix arg)
		{
			if(M!= arg.rows() || N != arg.cols())
			{
				MATRIX_OUT_ERR("Invalid operation mtx A = mtx B. Matrix dimension mismatch",__FILE__,__LINE__);
				throw std::invalid_argument("Matrix dimensions mismatch");
			}

			if(this != &arg)
				for(size_t i = 0; i<M; i++)
					for(size_t j = 0; j< N; j++)
						data[i*N + j ] = arg(i,j);


			return *this;
		}

		// make the function call operator retrive the i,j th data element
		complex double& operator()(int i,int j)
		{
			if(i >= M || i < 0 )
			{
				MATRIX_OUT_ERR("Array index out of bounds. Cannot retrieve cplxMatrix element. ",__FILE__,__LINE__)
					throw std::out_of_range("row index out of bounds.");
			}

			if(j >= N || j < 0 )
			{
				MATRIX_OUT_ERR("Array index out of bounds. Cannot retrieve cplxMatrix element. ",__FILE__,__LINE__)
					throw std::out_of_range("column index out of bounds.");
			}


			return  data[i*N+j];

		}

		complex double* data(){ return  data;}
		const size_t rows() const { return M; };
		const size_t cols() const { return N; };
		// overload the multiplication by matrix operator...
		// multiply an MxN by an NxK matrix and store the result in an MxK matrix C! 
		cplxMatrix operator*(cplxMatrix B)
		{

			if(B.rows() != N)
			{
				MATRIX_OUT_ERR("Matrix dimensions do not comply with product rule",__FILE__,__LINE__);
				throw std::invalid_argument("invalid argument exception!");
			}

			size_t K = B.cols();
			cplxMatrix C(M,K);

			complex double alpha = 333I*0;
			complex double beta = 0;
			cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M,K,N, &alpha,  data,M,B.data(),N, &beta,(void*) C.data(),M);
			return C;
		}
	};


	/**
	* return the identity matrix of size NxN
	*/
	mtx::cplxMatrix eye(size_t N){

		mtx::cplxMatrix res(N,N);
		for(size_t i = 0; i< N;i++)
			res(i,i) = 1. ;

		return res;
	}

	/**
		Return a ones matrix of size MxN
	*/
	mtx::cplxMatrix ones(size_t M,  size_t N){

		mtx::cplxMatrix res(M,N);
		for(size_t i = 0 ; i < M ; i++)
			for (size_t j = 0 ; j < N  ;j++)
				res(i,j)  = 1.;
		return res;
	}

}



#endif

