/* Class for sparse matrix diagonalization; crs used for storing the 
   matrix. */

#ifndef MATRIX_SPARSE_REAL_HPP
#define MATRIX_SPARSE_REAL_HPP


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
double GoToDot(int *,int *,double[],double[],double[],int *,int[],int *);
void GoToImage(int *,int *,double[],double[],double[],int *,int[],int *);
void GoToMonit(int *,int *,int *,int *,int *,double[],double[]);




class Matrix_sparse_real
{

  static const double DIAG_TOLERANCE; // precision at diagonalization
  static const int DIAG_MAX_ITERATIONS=1000;// max.Nr. of iterations at diagon.

  int *row_NrNz;     // Nr. of stored elements for each row
  int *col_ind;      // column indices of the stored elements
  double *mat_sp;    // value of the matEl
  int nz_elements,dimension;  // actual Nr. of non-zero elements, dimension
  int how_many_eigvects;      // ... to be found at diagonalization
  int add_to_diag_workspace_sparse; // k=m+add_to_diag...
  int foomax(int a,int b);    // just find the maximum of a,b
public:
  void set_k(int k);          // set a new value of add_to_diag_workspace_sparse
  /* Routine calculating A*v, see doc. for f02fjf_. It virtually contains 
     the information about the matrix. */
  void image(int *,double[],double[]);

  /* Routine calculating u*v, see doc. for f02fjf_ */  
  double dot(int *,double[],double[]);

  /* Routine which may be used to monitor progress of the diagonalization, 
     see doc. for f02fjf_ */  
  void monit(int *,int *,int *,int *,int *,double[],double[]);

  Matrix_sparse_real(int *,int *,double *,int,int,int,int);

  virtual ~Matrix_sparse_real();

  /* Guess what... parameters: 
     x       eigenvectors, stored one after another in the 1D array
             (size has to be at least (k+4)*dimension???).
     eigval  eigenvalues
  */
  int diagonalize(double *x,int novecs,double *eigval);
};

#endif

















