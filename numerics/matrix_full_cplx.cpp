
/*!
 \file matrix_full_cplx.cpp
   \brief Implementation of functions 
   for diagonalizing complex hermitean full matrices.
 */

#include <numerics/matrix_full_cplx.hpp>
#include <utils/logger.hpp>
#include <LLLlib/ERRORS.h>

#ifdef _NAG_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#else
#ifdef  INTEL
#include <mkl.h>
#include <mkl_lapack.h>
#endif
#ifdef GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#endif
#endif
//#include <atlas_enum.h>
//#include "clapack.h"

/*!
  \fn int diag_full_cplx_matrix(int dimension, int mtmp, dComplex *a, double *w, dComplex *z)
  \param dimension obvoius
  \param eigsToFind Number of requested eigenvalues
  \param *a 
  \param *w
  \param *z
  \brief diagonalizes a complex matrix, using the NAG library
  See the NAG-library documentation for explanation of 
   f02hcf_()-related symbols. 
*/
int diag_full_cplx_matrix(int dimension, int eigsToFind, dComplex *a,
			  double *w, dComplex *z)
{ 
#ifdef GSL
	#ifdef LAPACK
		return (diag_full_cplx_matrix_lapack(dimension,eigsToFind, a,w,z));
#else
	return (diag_matrix_gsl(dimension,eigsToFind,a,w,z));
#endif
	//return (diag_full_cplx_matrix_lapack(dimension,eigsToFind, a,w,z));
#else
#ifdef LAPACK
  return (diag_full_cplx_matrix_lapack(dimension,eigsToFind, a,w,z));
#else
  char 
  	job = 'V',
  	range ='I',
  	uplo = 'L';
  	
	
  //int blocksize = ilaenv(1,"ZHEEVR","VIL",&n, -1, -1,-1);
#ifdef MKL
  MKL_INT
    	n = dimension,
    	lda = dimension,
    	il = 1,
    	m = eigsToFind,
    	iu = eigsToFind,
    	mest = eigsToFind,
    	lwork = 2*dimension,
    	ifail = 0,
    	ldz=n,
    	isuppz = 2*eigsToFind,
    	lrwork = 70*n,
    	liwork = 50*n;
  double *rwork=new double[lrwork];
   MKL_INT * iwork=new MKL_INT[liwork];
   double wl =0.0,
     			wu = 0.0;

#else
int
  	n = dimension,
  	lda = dimension,
  	il = 1,
  	m = eigsToFind,
  	iu = eigsToFind,
  	mest = eigsToFind,
  	lwork = 2*dimension,
  	ifail = 0,
  	ldz=n,
  	isuppz = 2*eigsToFind,
  	lrwork = 70*n,
  	liwork = 50*n;

  double wl =0.0,
  			wu = 0.0;
  
 
  
  double *rwork=new double[lrwork];
  int * iwork=new int[liwork];
#endif

  dComplex *work=new dComplex[lwork];
  ifail=-1;
  /* How the matrix elements are stored in a...
     a[0].real=1;a[0].imag=0;            // a_11
     a[1].real=1.1;a[1].imag=1;          // a_21
     a[2].real=1.1;a[2].imag=-1;         // a_12
     a[3].real=2.1221;a[3].imag=0;       // a_22
  */
  #ifdef NAG
  f02hcf_(&job,&range,&uplo,&n,a,&lda,&wl,&wu,&il,&iu,&mest,&m,w,z,&ldz,work,&lwork,rwork,iwork,&ifail);
  #else 
	#ifdef MKL
	
	char prec ='S';
	double abstol = 2* DLAMCH(&prec);
	glLogger.error("Calling INTEL diag routine with (%f)", abstol);
	MKL_INT info = 0;
//	zheevx(&job,&range,&uplo,&n,a,&lda,&wl,&wu,&il,&iu,&abstol,&m,w,z,&ldz,work,&lwork,rwork,iwork,&ifail, &info);
	zheevr(&job,&range,&uplo,&n,a,&lda,&wl,&wu,&il,&iu,&abstol,&m,w,z,&ldz,
			&isuppz,work,&lwork,rwork,&lrwork, iwork, &liwork, &info);
	//! TODO handle info here!
	glLogger.error("... back from INTEL with (%d) and workspace suggestion (%d)", info);
	ifail = info;
	dComplex *temp = z;
	unsigned int count = 0;

	glLogger.debug("Leaving diagonalization with MKL");
	 delete[] work;
	  delete[] rwork;
	  delete[] iwork;
	  return 0;
          #else
	glLogger.error("Not implemented yet!! GSL");
    		throw (CxErrors(__FILE__,__LINE__));
		#endif
	#endif
#endif


#endif
	       
}

/*!
  \fn int calculateComplexDeterminant (int dimension, dComplex *theMatrix, dComplex *deterMinante)
*/
int calculateComplexDeterminant (int dimension, 
				 dComplex *theMatrix, 
				 std::complex<double> & determinante)
{
  std::complex<double> nullElement(0.0, 0.0);
  if (theMatrix == 0)
    {
      determinante = nullElement;
      return (-1);
    }
    if (dimension == 2)
    {
    	int retVal = calculateTwoDeterminante(theMatrix,determinante);
    	if (retVal != 0)
    	{
    		throw (CxErrors(__FILE__,__LINE__));
    	}
    	return (retVal);
    }
    Eigen::MatrixXcd daMatrix(dimension,dimension);
    for (int rowIndex = 0; rowIndex < dimension; ++rowIndex)
    		{
    			for (int colIndex = 0; colIndex < dimension; ++colIndex)
    			{
    				/*
    				gsl_complex value = gsl_complex_rect(theMatrix->real, theMatrix->imag);
    				gsl_matrix_complex_set(matrix, colIndex, rowIndex, value);
    				theMatrix++;
    				*/
    				std::complex<double> tempVal;
    				tempVal.real(theMatrix->real());
    				tempVal.imag(theMatrix->imag());
    				daMatrix(rowIndex,colIndex) = tempVal;

    				theMatrix++;
    			}
    		}
    determinante = daMatrix.determinant();
    return 0;
#ifdef NAG
  int 
    iFail,
    myDim;

  myDim = dimension;
  iFail = -1;
  

double 
    realDet = 0.0,
    imagDet = 0.0;
 double *workSpace = new double [dimension*2];
  f03adf_(theMatrix,
	  &myDim,
	  &myDim,
	  &realDet,
	  &imagDet,
	  workSpace,
	  &iFail);
	   delete [] workSpace;  
  switch (iFail)
    {
    case (1): 
      determinante = std::complex<double> ( realDet,imagDet);
      glLogger.debug("calculateComplexDeterminate With NAG library return (%f) + i(%f)", realDet, imagDet);
      return (0);
      break;
    case (2):
      return (2);
      break;
    case (3):
      return (3);
      break;
    case (4):
      return (4);
      break;
    }
  determinante = std::complex<double> ( realDet,imagDet);
	  #endif 
	 #ifdef GSL
		return (calculateComplexDeterminantGsl(dimension,theMatrix,determinante));
	
	  #endif
	 #ifdef MKL
	  return (calculateComplexDeterminanteMKL(dimension,theMatrix,determinante));
	 #endif
	
 
  //  delete realDet

  return 0;
}


/*!
  \fn int calculateComplexDeterminant (int dimension, dComplex *theMatrix, dComplex *deterMinante)
*/
#ifdef GSL
int calculateComplexDeterminantGsl (int dimension, 
				 dComplex *theMatrix, 
				 std::complex<double> & determinante)
{
  std::complex<double> nullElement(0.0, 0.0);
  if (theMatrix == 0)
    {
      determinante = nullElement;
      return (-1);
    }
    if (2 ==dimension)
    {
    	return calculateTwoDeterminante(theMatrix, determinante);
    }
   
   gsl_matrix_complex * m_matrix = gsl_matrix_complex_calloc (dimension, dimension);
   allocateAndFillMatrix(m_matrix,dimension,theMatrix);
   // fill matrix with values from theMatrix
	
// LU decomposition needed before determinante n berechnung
	gsl_permutation *m_perm = gsl_permutation_alloc(dimension);
	int signum = 0;
	
	gsl_linalg_complex_LU_decomp(m_matrix, m_perm, &signum);
	// calculate determinante 
	
	  gsl_complex m_determinante = gsl_linalg_complex_LU_det(m_matrix, signum);
	  determinante = std::complex<double>(GSL_REAL(m_determinante), GSL_IMAG(m_determinante));
	    glLogger.debug("calculateComplexDeterminateGsl return (%f) + i(%f)", determinante.real(), determinante.imag());
	  gsl_matrix_complex_free(m_matrix);
	  gsl_permutation_free(m_perm);
 
  //  delete realDet

  return 0;
}


#endif
int calculateTwoDeterminante(dComplex *theMatrix, std::complex <double> & determinante)
{
	#ifdef GSL
	gsl_complex val0=gsl_complex_rect(theMatrix->real, theMatrix->imag);
	theMatrix++;
	gsl_complex val1=gsl_complex_rect(theMatrix->real, theMatrix->imag);
	theMatrix++;
	gsl_complex val2=gsl_complex_rect(theMatrix->real, theMatrix->imag);
	theMatrix++;
	gsl_complex val3=gsl_complex_rect(theMatrix->real, theMatrix->imag);
	
	
	gsl_complex tempDet = gsl_complex_sub(gsl_complex_mul(val0,val2) , gsl_complex_mul(val1,val3));
	glLogger.debug("calculateTwoDeterminante returns (%f),(%f)",GSL_REAL(tempDet),GSL_IMAG(tempDet));
	determinante = std::complex<double>(GSL_REAL(tempDet),GSL_IMAG(tempDet));
	return (0);
	#else
	// for MKL here
	std::complex<double> val0(theMatrix->real(), theMatrix->imag());
	theMatrix++;
	std::complex<double> val1(theMatrix->real(), theMatrix->imag());
	theMatrix++;
	std::complex<double> val2(theMatrix->real(), theMatrix->imag());
	theMatrix++;
	std::complex<double> val3(theMatrix->real(), theMatrix->imag());
	std::complex<double> tempDet = (val0*val2)-(val1*val3) ;
	glLogger.debug("calculateTwoDeterminante returns (%f),(%f)",tempDet.real(),tempDet.imag());
	determinante = tempDet;
	return (0);	
	#endif
}
							 	
							 	
							 
int calculateThreeDeterminante(dComplex *theMatrix, std::complex <double> & deterMinante)
{
	std::complex<double> retVal;
	
	return (-1);	
}

int calculateComplexDeterminanteMKL(int dimension,dComplex *theMatrix, std::complex <double> &determinante)
{
	#ifdef MKL
	// Do the same as in GSL: LU factorization and readout of determinante
	 MKL_INT mklDimension = dimension;
	MKL_INT *ipiv = new MKL_INT[dimension];
	MKL_INT info = 0;
/* From Include file	
 * void zgetrf_(int *m,int *n, MKL_Complex16 *a,int *lda,int *ipiv,int *info);
 */
 if (glLogger.getLogLevel() == DEBUG)
 {
 	glLogger.info("Before LAPACK call");
 	dComplex *tempMatPtr = theMatrix;
 	for (int colIndex = 0; colIndex < dimension; ++colIndex) 
 	{
		for (int rowIndex = 0; rowIndex < dimension; ++rowIndex) 
		{
			// show matrix
 			glLogger.debug("Matrix: element (%d), (%d) = (%f) + i(%f)", rowIndex, colIndex, tempMatPtr->real(),tempMatPtr->imag());
 			tempMatPtr++;	
		}
	}
 	
 }

 
 
	 zgetrf(&mklDimension, &mklDimension, theMatrix,&mklDimension, ipiv, &info);
	// theMatrix is changed now, we read the diagonal elements
 if (glLogger.getLogLevel() == DEBUG)
 {
 	glLogger.debug("After LAPACK call");
 	dComplex *tempMatPtr = theMatrix;
 	for (int colIndex = 0; colIndex < dimension; ++colIndex) 
 	{
		for (int rowIndex = 0; rowIndex < dimension; ++rowIndex) 
		{
			// show matrix
 			//glLogger.debug("Matrix: element (%d), (%d) = (%f) + i(%f)", rowIndex, colIndex, tempMatPtr->real,tempMatPtr->imag);
 			tempMatPtr++;	
		}
	}
 	
 }


	std::complex<double> storage(1.0,0.0);
		for (unsigned int pointer = 0; pointer < dimension*dimension; pointer=pointer + dimension +1) 
	{
		
		float tempReal = theMatrix[pointer].real();
		float tempImage = theMatrix[pointer].imag();
		std::complex<double> tempElement(tempReal,tempImage);
		storage = storage * tempElement;
	}
	int sign = -1;
	for (int count = 0; count < dimension; ++count) {
		glLogger.debug("Ipiv Matrix element (%d) = (%d)", count, ipiv[count]);
		if (ipiv[count] != count)
		{
			// Trial here, see LAPACK docu for more info
			sign = sign *-1;
		}
	}
	determinante = storage*std::complex<double>(sign,0.0);
	
	glLogger.info("Factorization from calculateComplexDeterminanteMKL returns (%d), determinante is (%f)+i(%f)", info, determinante.real(), determinante.imag());
	return (0);
	// read out determinante 
	#else
		ERRORTHROW("Not implemented");
	#endif
	return (-1);	
}
#ifdef MKL
int multiplyComplex(dComplex& factor1, dComplex &factor2, dComplex & factor)
{
	// factor1 = a+ib, factor 2 = c+id
	double a = factor1.real();
	double b = factor1.imag();
	double c = factor2.real();
	double d = factor2.imag();
	double facReal = (a*c) - (b*d);
	double facImag = (a*d) + (b*c);
	factor.real(facReal);
	factor.imag (facImag);
return (0);	
}
#endif


int diag_matrix_gsl(int dimension,int eigsToFind,dComplex *a,double *w,dComplex *z)
{
#ifndef GSL
	glLogger.error("GSL diagonalization called but not implemented!");
	throw (CxErrors(__FILE__,__LINE__));

#else 
	glLogger.error("GSL diagonalization called with dimension %d", dimension);
	gsl_matrix_complex *matrix = gsl_matrix_complex_calloc (dimension, dimension);
	if (matrix == 0)
	{
		
		THROWERROR;
	}
	if (!allocateAndFillMatrix(matrix,dimension,a))
	{
		THROWERROR;
		
	}
	
	    
       gsl_vector *eval = gsl_vector_alloc (dimension);
       gsl_matrix_complex *evec = gsl_matrix_complex_alloc (dimension, dimension);
     
       gsl_eigen_hermv_workspace * workspace = 
         gsl_eigen_hermv_alloc (dimension);
       
       gsl_eigen_hermv (matrix, eval, evec, workspace);
       ERROR("START Now with diagonalization");
       gsl_eigen_hermv_free (workspace);
     
       gsl_eigen_hermv_sort (eval, evec,GSL_EIGEN_SORT_VAL_ASC);

       gsl_matrix_complex_free(matrix);
       
       
       dComplex *arrayPointer = z;
       // Lenght of a line = dimension * sizeof dcomplex
     
       
       // copy eigenvalues and eigenstates
       
       for (int evalIndex = 0; evalIndex < eigsToFind;evalIndex++)
       {
    	// get eigenvalue
    	   double eval_i  = gsl_vector_get (eval, evalIndex);
    	   w[evalIndex] = eval_i;
    	 
    	   glLogger.error ("eigenvalue No (%d) = %g\n",evalIndex, eval_i);
    	 //  glLogger.error ("eigenvector = \n");
    	if (evalIndex == 0) 
    	{
    		gsl_vector_complex_view evec_i =  gsl_matrix_complex_column (evec, evalIndex);

    	   gsl_vector_complex_fprintf (stdout, &evec_i.vector, "%g");
    	}  
    	
    	   for (int evecIndex = 0;evecIndex < dimension;  evecIndex++)
    	   {
    		   gsl_complex temp = gsl_matrix_complex_get(evec,evecIndex, evalIndex);
    		   
    		   arrayPointer->real = GSL_REAL(temp);
    		   arrayPointer->imag = GSL_IMAG(temp);
    		   
    		   //glLogger.debug("evec[%d,%d] = (%f)=(%f), (was (%f , %f)",evalIndex,evecIndex,arrayPointer->real, arrayPointer->imag, GSL_REAL(temp), GSL_IMAG(temp));
    		   arrayPointer ++;
    		// get coefficents   
    	   }
    	   
    	   
    	   
    	   
       }
       
       
       gsl_matrix_complex_free(evec);
      
       
       
	return (dimension);











#endif


}

int diag_full_cplx_matrix_lapack(int dimension,int eigsToFind,dComplex *a,double *w,dComplex *z)
{
#ifdef LAPACK
	ERROR("Diagonalization with LAPACK");
	char 
	  	job = 'V',
	  	range ='I',
	  	uplo = 'L';
	  	
		int n = dimension;
/*	dont understand what ilaenv is supposed to do may come back later	
int 
			n1 =1,
			n2 = -1,
			n3 = -1,
			n4 = -1;
	  int blocksize = ilaenv_(&n1,"ZHEEVR","VIL",&n, &n2, &n3,&n4);
*/
	int
	   	lda = dimension,
	  	il = 1,
	  	m = eigsToFind,
	  	iu = eigsToFind,
	  	lwork = 2*dimension,
	  	ifail = 0,
	  	ldz=n,
	  	isuppz = 2*eigsToFind,
	  	lrwork = 70*n,
	  	liwork = 50*n;

	  double wl =0.0,
	  			wu = 0.0;
	  
	 
	  
	  double *rwork=new double[lrwork];
	  int * iwork=new int[liwork];

	  dComplex *work=new dComplex[lwork];
	  ifail=-1;
	  /* How the matrix elements are stored in a...
	     a[0].real=1;a[0].imag=0;            // a_11
	     a[1].real=1.1;a[1].imag=1;          // a_21
	     a[2].real=1.1;a[2].imag=-1;         // a_12
	     a[3].real=2.1221;a[3].imag=0;       // a_22
	  */
	  		
		char prec ='S';
		double abstol = 2* dlamch_(&prec);
		glLogger.error("Calling LAPACK diag routine with (%f)", abstol);
		int info = 0;
	
		zheevr_(&job,&range,&uplo,&n,a,&lda,&wl,&wu,&il,&iu,&abstol,&m,w,z,&ldz,
				&isuppz,work,&lwork,rwork,&lrwork, iwork, &liwork, &info);
		//! TODO handle info here!
		glLogger.error("... back from LAPACK with (%d) and workspace suggestion (%d), ifail = (%d)", info, work[0],ifail);
		ifail = info;

		
	
return (1);	
	
#else 
THROWERROR("SHOULD NOT BE HERE");
return (-1);

#endif

	
	
}



#ifdef GSL

bool allocateAndFillMatrix(gsl_matrix_complex * matrix, int dimension, dComplex *theMatrix)
{
	
     //matrix = gsl_matrix_complex_calloc (dimension, dimension);
	// fill matrix with values from theMatrix
		for (int rowIndex = 0; rowIndex < dimension; ++rowIndex) 
		{
			for (int colIndex = 0; colIndex < dimension; ++colIndex) 
			{
				gsl_complex value = gsl_complex_rect(theMatrix->real, theMatrix->imag);
				gsl_matrix_complex_set(matrix, colIndex, rowIndex, value);
				theMatrix++;
			}
		}
	
	return (true);
}


#endif
