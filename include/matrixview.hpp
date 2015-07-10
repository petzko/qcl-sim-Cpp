/*
 * matrixview.hpp
 *
 *  Created on: Jul 10, 2015
 *      Author: petzko
 */

#ifndef INCLUDE_MATRIXVIEW_HPP_
#define INCLUDE_MATRIXVIEW_HPP_
#include <matrix.hpp>



/**
 * The matrix view class shall be a subclass of Matrix that creates a view window of a subset of the parent matrix's elements.
 *
 * Important characteristics of matrix view:
 *
 * 1) Each operation possible on a matrix shall be possible on a view! However the results of this operations SHALL be updated in the parent matrix as well!!
 * 2) the view shall not own the matrix data.
 * 3) the view shall be in general shorter lived than the parent matrix. this means that whenever view is destroyed the parent matrix shall NOT be affected,
 * however whenever a parent matrix is destroyed, all of its corresponding views shall be destroyed.
 *
 * The implementation will be based on the observer pattern! Each view subscribes to the parent and gets assigned a unique subscription id. This id shall be used later upon a view object destruction
 * to notify the parent object of which view exactly is destroyed. Destruction of the view shall leave the state of the parent object UNALTERED!!!!
 *
 *	ToDo: Test other options and decide which one is faster
 *	implementation strategy:
 *  ---- as all operations on the view shall effect the state of the parent one can implement those operations (still using the BLAS routines) onto a copy of the parent data followed by an inverse update
 *  of the parent data with the corresponding changes made to the view! Note that the general Level 1 and Level 2 BLAS routines allow for noncontiguous memory operations -> for example the routine
 *  xAXPY, which we use for vector (matrix addition) has the following signature:
 *  xAXPY(N,ALPHA,X,INCX,Y,INCY) -> where INCX and INCY can be chosen accordingly when operating on a non-contiguous chunk of
 *
 *  ---- For example. Given the parent NxN matrix M and the view view = M(1,N-1,1,N-1); and the scaling operation view = view*2; internally this scaling will be implemented as follows:
 *  -----
 */

namespace dat{
template<typename _Tp>
class MatrixView: public Matrix<_Tp> {

private:

	Matrix<_Tp>* _parent;
	unsigned int _id;

	unsigned int _dimI, _dimJ, _bytes;
	_Tp** _viewdata;

public:

	MatrixView(Matrix<_Tp>* parent, unsigned int id, unsigned int i1, unsigned int i2, unsigned int j1, unsigned int j2);

	~MatrixView();
//
//	void operator=(const Matrix<_Tp>& arg);
//
//	/**
//	 * Return a reference (implicit pointer) to the i,j-th data element
//	 * so that the user can actually overwrite the corresponding entry
//	 */
//	_Tp& operator()(int i, int j) const;
//
//	/**
//	 * return an error - cannot create a view of the view!
//	 */
//
//	Matrix<_Tp> operator()(int i1, int i2, int j1, int j2);
//	// overload the multiplication by matrix operator...
//	Matrix<_Tp> operator*(const Matrix<_Tp>& arg) const;
//	//overload the addition operator
//	Matrix<_Tp> operator+(const Matrix<_Tp>& arg) const;
//	//overload the subtraction operator
//	Matrix<_Tp> operator-(const Matrix<_Tp>& arg);
};

}


#endif /* INCLUDE_MATRIXVIEW_HPP_ */
