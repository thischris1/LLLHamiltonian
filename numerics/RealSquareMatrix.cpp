/*! 
 \file RealSquareMatrix.cpp
 \brief Implementation of file
 */
 
#include <numerics/RealSquareMatrix.h>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>


RealSquareMatrix::~RealSquareMatrix()
{
	
}

RealSquareMatrix::RealSquareMatrix(int rank):
SquareMatrix(rank),
m_matrix(0)
{
	m_matrix = gsl_matrix_calloc(rank, rank);
}

RealSquareMatrix::RealSquareMatrix(int rank, double *n_data):
SquareMatrix(rank),
m_matrix(0)
{
	m_matrix = gsl_matrix_calloc(rank, rank);
	int count = 0;
	for (int rowIndex = 0; rowIndex < rank; ++rowIndex) {
		for (int colIndex = 0; colIndex < rank; ++colIndex) {
			gsl_matrix_set(m_matrix,rowIndex, colIndex, n_data[count]);
			count++;
		}
		
	}
		
}

CMDCOMPLEX RealSquareMatrix::getElement(const int row, const int col)
{
	CMDCOMPLEX retVal(0.0,0.0);
	(gsl_matrix_get(m_matrix, row, col), 0.0);
	return (retVal);	
	
}

bool RealSquareMatrix::setElement(const CMDCOMPLEX n_value, const int row, const int col)
{
	double n_val = n_value.real();
	gsl_matrix_set(m_matrix,row, col, n_val);
	gsl_matrix_set(m_matrix,col, row, n_val);
	return (true);
	
}
int RealSquareMatrix::getSpectrum(DOUBLE_VECTOR &eigVals, COMPLEX_VECTOR &eigVecs)
{
	eigVals = DOUBLE_VECTOR(getCols());

	gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (getCols());
	gsl_vector *eval = gsl_vector_alloc (getCols());
    gsl_matrix *evec = gsl_matrix_alloc (getCols(), getCols());
	gsl_eigen_symmv (m_matrix, eval, evec, workspace );
	gsl_eigen_symmv_free(workspace);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	/* Copy results
	 * */

	for (int index = 0; index < getCols(); ++index) 
	{
		eigVals[index]=gsl_vector_get(eval, index);
		gsl_vector_view evec_i = gsl_matrix_column (evec, index);
		for (int vecIndex = 0; vecIndex < getCols(); ++vecIndex) {
			
			MDOUBLE dval = gsl_vector_get(&evec_i.vector, vecIndex);
			eigVecs[index*getCols()+vecIndex]=CMDCOMPLEX(dval,0.0);
		}
				
	}

	return (0);	
	
}

