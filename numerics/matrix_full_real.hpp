/* Function for diagonalizing real symmetric full matrices.
 */

#ifndef MATRIX_FULL_REAL_HPP
#define MATRIX_FULL_REAL_HPP 

// NAG fortran library routine for diagonalization of a real sym. matrix
#ifdef NAG
extern "C" void f02fcf_(char *,char *,char *,int *,double[],int *,double *,double *,int *,int *,int *,int *,double[],double[],int *,double[],int*,int[],int *);
#endif
/* Arguments:
   dimension  self-explaining
   m          how many (lowest?) eigenvalues+eigvec to be found
   a          the matrix (stored row-wise as doubles, i.e.
              a_00,a_01,a_02,...); only the upper triangle has to be set
   w          for eigenvalues
   z          for eigenvectors (stored the same as el's in a)
*/
int diag_full_real_matrix(int dimension,int m,double *a,double *w,double *z);


#endif
