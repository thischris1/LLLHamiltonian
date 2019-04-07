#ifndef REALSQUAREMATRIX_H_
#define REALSQUAREMATRIX_H_

#include "SquareMatrix.h"
#include <LLLlib/LLLlib.h>
#include <gsl/gsl_matrix.h>
class RealSquareMatrix : public SquareMatrix
{
public:
	RealSquareMatrix(int rank);
	RealSquareMatrix(int rank, double* value);
	
	virtual ~RealSquareMatrix();
	CMDCOMPLEX getElement(const int row, const int col);
	int getSpectrum(DOUBLE_VECTOR & eigVals, COMPLEX_VECTOR & eigVecs);
	bool setElement(const CMDCOMPLEX  n_value, const int row, const int col);
private:
	// helper
	// do not implement
	RealSquareMatrix();
	// members

	gsl_matrix *m_matrix;

};

#endif /*REALSQUAREMATRIX_H_*/
