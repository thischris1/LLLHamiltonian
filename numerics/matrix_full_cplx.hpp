/*! 
  \file  matrix_full_cplx.hpp
   \brief Declaration of Function for diagonalizing complex hermitean full matrices.
 */

#ifndef MATRIX_FULL_CPLX_HPP
#define MATRIX_FULL_CPLX_HPP 
#include <complex>
#include <LLLlib/Types.hpp>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#ifdef NAG
//! NAG fortran library routine for diagonalization of a compl. herm. matrix
extern "C" void f02hcf_(char *,char *,char *,int *,dComplex[],int *,double *,double *,int *,int *,int *,int *,double[],dComplex[],int *,dComplex[],int*,double[],int[],int *);

//! NAG Fortran library routine for calculation of determinant of a complex matrix. See NAG docu for details
extern "C" void f03adf_(dComplex[] , int *, int *, double *, double *, double *, int *);

//! NAG Fortran library routine for minimizing a function

#endif




int findMinimaforComplexFunction();

/*! 
  \fn int diag_full_cplx_matrix(int dimension,int m,dComplex *a,double *w,dComplex *z);
  \param  dimension  self-explaining
  \param  m  how many (lowest?) eigenvalues+eigvec to be found
  \param  a the matrix (stored row-wise as dComplex'es, i.e.
  a.re_00,a.im_00,a.re_01,a.im_01,a.re_02,...
  only the upper triangle has to be set
  \param w          for eigenvalues
  \param z          for eigenvectors (stored the same as el's in a)
*/

int diag_full_cplx_matrix(int dimension,int m,dComplex *a,double *w,dComplex *z);
int diag_full_cplx_matrix_lapack(int dimension,int m,dComplex *a,double *w,dComplex *z);
int diag_matrix_gsl(int dimension,int m,dComplex *a,double *w,dComplex *z);

/*! 
  \fn int calculateComplexDeterminant (int dimension, dComplex *theMatrix, dComplex *deterMinante)
  \param *deterMinante  the determinante (on return)
  \param dimension Guess what
  \param *theMatrix the Matrix to be determined
  \return the determinante in Determinante, the errorcode of NAG aus return Value
  
  
*/
int calculateComplexDeterminant (int dimension, 
				 dComplex *theMatrix, 
				 std::complex <double> & deterMinante);
		      
int calculateComplexDeterminantGsl (int dimension, 
							 dComplex *theMatrix, 
							 std::complex <double> & deterMinante);

int calculateTwoDeterminante(dComplex *theMatrix, 
							 std::complex <double> & deterMinante);
							 
int calculateThreeDeterminante(dComplex *theMatrix, 
							 std::complex <double> & deterMinante);
							 
int calculateComplexDeterminanteMKL(int dimension, 
							 dComplex *theMatrix, 
							 std::complex <double> & );


#ifdef MKL
// convenience routine
	int multiplyComplex(dComplex& factor1, dComplex &factor2, dComplex & factor);
#endif
#ifdef GSL
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_complex.h>
	bool allocateAndFillMatrix(gsl_matrix_complex * matrix, int dimension, dComplex *theMatrix);
#endif 

#ifdef LAPACK
	extern "C" void    zheevr_( char *jobz, char *range, char *uplo, int  *n, dComplex *a, int  *lda, double *vl, double *vu, int  *il, int  *iu, double *abstol, int  *m, double *w, dComplex *z, int  *ldz, int *isuppz, dComplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int * liwork, int * info );
	
	extern "C" double dlamch_(char * s);
	extern "C" int  ilaenv_(  int *ispec, char *name,  char *opts,  int *n1,  int *n2, int *n3, int *n4 );
#endif
#endif
