/* Class for sparse matrix diagonalization; crs used for storing the 
   matrix. */

#ifdef COMPILE_PARALLELIZED

#ifndef MATRIX_SPARSE_REAL_PARALLEL_HPP
#define MATRIX_SPARSE_REAL_PARALLEL_HPP




#include<vector>
#include<mpi_tag_code.hpp> 
#include<mpi.h>
#include<logger.hpp>
#include<assert.h>


// NAG fortran library routine for diagonalization of a sparse real sym. matrix
extern "C" void f02fjf_(int *,int *,int *,int *,double *,double (*dot)(int *,int *,double[],double[],double[],int *,int[],int *),void (*image)(int *,int *,double[],double[],double[],int *,int[],int *),void (*monit)(int *,int *,int *,int *,int *,double[],double[]),int *,double[],int *,double[], double[],int *,double[],int *,int[],int *,int *);

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
double GoToDot_parallel(int *,int *,double[],double[],
			double[],int *,int[],int *);
void GoToImage_parallel(int *,int *,double[],double[],
			double[],int *,int[],int *);
void GoToMonit_parallel(int *,int *,int *,int *,
			int *,double[],double[]);




class Matrix_sparse_real_parallel
{
public:
  // structure for sparse storing of matrix els.
  typedef struct {
    double val;
    int col;
  } MatEl;

private:
  static const double DIAG_TOLERANCE; // precision at diagonalization
  static const int DIAG_MAX_ITERATIONS=10000;// max.Nr. of iterations at diagon.

  int *p_row_NrNz;     // Array of Nr. of stored elements for each row
  std::vector<MatEl> **p_rowOfMatEls;  // column indices + values of the stored elements
  int nz_elements;   // actual Nr. of non-zero elements (within the current proc)
  int dimension;     // dimension 
  int rowEvaluated;  // Nr. of rows handled by the current proc
  int how_many_eigvects;      // ... to be found at diagonalization
  int add_to_diag_workspace_sparse; // k=m+add_to_diag...
  int eigsToFind;

  int numprocs;
  int myid;
  //MPI_Request *req_productReady;
  double *vecTmp;    // for storing received vectors

  int foomax(int a,int b);    // just find the maximum of a,b
public:
  
  void set_k(int k);          // set a new value of add_to_diag_workspace_sparse
  /* Calculation of A*v for the master proc. Collects results of the slave
     procs calculated in image_part. */
  void image(int *iflag,double z[],double w[]);

  /* Routine calculating a contribution to A*v, see doc. for f02fjf_. The whole
     A*v is the sum of outVec from all slave processes. This routine 
     virtually contains the information about the matrix. */
  void image_part(double inVec[],double outVec[]);

  /* Routine calculating u*v, see doc. for f02fjf_ */  
  double dot(int *,double[],double[]);

  /* Routine which may be used to monitor progress of the diagonalization, 
     see doc. for f02fjf_ */  
  void monit(int *,int *,int *,int *,int *,double[],double[]);

  Matrix_sparse_real_parallel(int *new_p_row_NrNz,
			      std::vector<MatEl> **new_p_rowOfMatEls,
			      int new_nz_elements,
			      int new_dimension,
			      int new_totRowEvaluated,
			      int new_eigsToFind,
			      int new_add_to_diag_workspace);

  virtual ~Matrix_sparse_real_parallel();

  /* Guess what... parameters: 
     x       eigenvectors, stored one after another in the 1D array
             (size has to be at least (k+4)*dimension???).
     eigval  eigenvalues
  */
  int diagonalize(double *x,int novecs,double *eigval);
};

#endif

#endif















