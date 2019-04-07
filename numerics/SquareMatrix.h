#ifndef SQUAREMATRIX_H_
#define SQUAREMATRIX_H_

#include "Matrix.h"


class SquareMatrix : public Matrix
{
public:
	SquareMatrix( int rank);
	virtual ~SquareMatrix();
//! Diagonalization takes place here
	virtual int getSpectrum(DOUBLE_VECTOR & eigVals, COMPLEX_VECTOR &eigVecs)=0;
	
private:
	SquareMatrix();
};

#endif /*SQUAREMATRIX_H_*/
