/*!
  \file Matrix_sparse_cplx.hpp
  \author Christian Mueller
  \brief declaraation of class Matrix_sparse_cplx
*/ 

#ifndef MATRIX_SPARSE_CPLX_HPP
#define MATRIX_SPARSE_CPLX_HPP
#include <complex>
#include <vector>
#include <fstream>
#ifdef MKL_ILP64
#include "mkl.h"
#endif
// NAG fortran library routine for diagonalization of a sparse real sym. matrix
#ifdef NAG
extern "C" void f02fjf_(int *,int *,int *,int *,double *,double (*dot)(int *,int *,double[],double[],double[],int *,int[],int *),void (*image)(int *,int *,double[],double[],double[],int *,int[],int *),void (*monit)(int *,int *,int *,int *,int *,double[],double[]),int *,double[],int *,double[], double[],int *,double[],int *,int[],int *,int *);
#endif
/* Mirrors of dot,image and monit (see below). Purpose: pointers to 
functions which are to be passed to f02fjf_() have to be of type 'pointer 
to a function' and not 'pointer to a method of class Matrix_sparse_real'. 
A hidden extra argument would be added to the function arguments in the 
latter case; it would be 'this' (pointer to the object of the 
Matrix_sparse_real class the function belongs to).
   Additionally, since (at least) GoToImage() has to call really a member 
method of the class, image(), the value of 'this' is passed to GoToDot as 
int *iwork (second last call parameter of GoToImage()). Note that iwork 
thus does NOT point to an (array of) integers but to a single object of 
type Matrix_sparse_real. See the implementation of diagonalize() 
and GoToImage().
*/
double VecDotVec(int *,int *,double[],double[],double[],int *,int[],int *);
void MatrixDotVectorWrapper(int *,int *,double[],double[],double[],int *,int[],int *);
void GoToMonit(int *,int *,int *,int *,int *,double[],double[]);


/*! 
  \class Matrix_sparse_cplx
 \brief Class for complex sparse matrix diagonalization crs used for storing the  matrix. 
 Class for complex sparse matrix diagonalization crs used for storing the  matrix. 
*/
//! Struct to model a sparse matrix element
typedef struct {std::complex<double> matEl; int rowIndex; int colIndex;} sparseMatEl;
class Matrix_sparse_cplx
{
private:
  static const double DIAG_TOLERANCE; // precision at diagonalization
  static const int DIAG_MAX_ITERATIONS;// max.Nr. of iterations at diagon.
  static const int ADD_TO_DIAG_WORKSPACE_SPARSE;
  std::vector <MKL_INT> row_NrNz;     //! Nr. of stored elements for each row
  std::vector <MKL_INT> col_ind;      //! column indices of the stored elements
  std::vector <std::complex<double> >  mat_sp;    //! value of the matEl
  int nz_elements; //! actual Nr. of non-zero elements
    int dimension;  //!  dimension
  int how_many_eigvects;      //! ... to be found at diagonalization
  double m_m0;
  int add_to_diag_workspace_sparse; //! k=m+add_to_diag...
  int foomax(int a,int b);    //! just find the maximum of a,b
  bool calculateReal;
  bool useRCI;
  double offSetDiagonal; //! Since the used method calculates the largest abs. of  eigenvalues and eigenvectors but we are interested in the lowest part of the spectrum, we need to substract a thumbvalue to move the whole spectrum to the negative.
  double emin;
  double emax;
  //! Private default ctor
  Matrix_sparse_cplx();
  Matrix_sparse_cplx ( const Matrix_sparse_cplx &);
  Matrix_sparse_cplx & operator = (const Matrix_sparse_cplx &);

public:
  void set_k(int k);          // set a new value of add_to_diag_workspace_sparse
  /* Routine calculating A*v, see doc. for f02fjf_. It virtually contains 
     the information about the matrix. */
  void matrixDotVector(int *, double *z,double * w);

  //! Routine calculating u*v, see doc. for f02fjf_ */  
  double dot(int *,double[],double[]);

  /*! Routine which may be used to monitor progress of the diagonalization, 
     see doc. for f02fjf_ */  
  void monit(int *,int *,int *,int *,int *,double[],double[]);


  Matrix_sparse_cplx(std::vector<MKL_INT> ,std::vector<MKL_INT> ,std::vector<std::complex<double> >,int,int,int,int);

  virtual ~Matrix_sparse_cplx(); //! Destructor

  /*! Guess what... parameters: 
     \param *eigVec       eigenvectors, stored one after another in the 1D array
             (size has to be at least (k+4)*dimension???).
     \param *eigval  eigenvalues
  */
  int diagonalize(std::complex<double> *eigVec,
		  int novecs,
		  double *eigval);

  int diagLowestOnly(std::complex<double> *eigVec,
		     int novecs,
		     double *eigval);
  //! Write matrix to screen
  void dumpMatrix(void) const;

  void dumpMatrix(std::ostream &outStream) const;
  void dumpCSRMatrix(std::ostream &outStream) const;
private:
 
};

#endif

















