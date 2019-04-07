#ifndef MATRIX_SPARSE_COMPLEX_PARALLEL_H_
#define MATRIX_SPARSE_COMPLEX_PARALLEL_H_
#include <LLLlib/LLLlib.h>

#include <LLLlib/Types.hpp>
int calculateComplexDeterminant (int dimension, 
				dComplex *theMatrix, 
				 std::complex <double> & deterMinante);

class Matrix_sparse_complex_parallel
{
public:
	Matrix_sparse_complex_parallel();
	virtual ~Matrix_sparse_complex_parallel();
};

#endif /*MATRIX_SPARSE_COMPLEX_PARALLEL_H_*/
