 /* Function for diagonalizing real symmetric full matrices.
 */


#include<numerics/matrix_full_real.hpp>
#include <iostream>
#ifndef NAG
#ifdef MKL
#include "mkl.h"
#endif
#endif
// RUNTIME includes
#include <utils/logger.hpp>
/* See the NAG-library documentation for explanation of 
   f02fcf_()-related symbols. */
int diag_full_real_matrix(int dimension,int mtmp,double *a,double *w,double *z)
{
#ifndef LAPACK
  char job,range,uplo;
  int n,lda,il,iu,mest,lwork,ifail,ldz,m;
  int *iwork;
  double wl,wu;
  double *work;
  
  job='V';
  range='I';
  uplo='L';
  n=dimension;
  lda=n;
  m=mtmp;
  il=1;
  iu=m;
  mest=m;
  ldz=n;
  lwork=8*n;
  wu=2;wl=1;
  iwork=new int[5*n];
  work=new double[lwork];
  ifail=-1;
  int info = 0;
  double absTol = 1e-06;
#endif
  #ifdef NAG

  
  /*  How the matrix elements are stored in a...
      a[0]=1;            // a_11
      a[1]=1.1;          // a_21
      a[2]=1.1;          // a_12
      a[3]=2.1221;       // a_22
  */
  f02fcf_(&job,&range,&uplo,&n,a,&lda,&wl,&wu,&il,&iu,&mest,&m,w,z,&ldz,work,&lwork,iwork,&ifail);
  #else
#ifdef MKL
  int *isuppz = 0; // not needed according to spec
  glLogger.error("About to enter dsyevr in dsyevr");
  /*
  dsyevr(&job,&range,&uplo,(MKL_INT *)(&n),a,(MKL_INT *)(&lda),&wl,&wu,(MKL_INT *)(&il),
		  (MKL_INT *)(&iu),&absTol,(MKL_INT *)(&m),w,z,(MKL_INT *)(&ldz),(MKL_INT *)(isuppz),work,(MKL_INT *)(&lwork),iwork,
  &ifail, &info);*/
  glLogger.error("Back with (%d),(%d)", ifail, info);
#endif
#ifndef LAPACK

  

  std::cout << "Exit code: "<< ifail << "\n";
  /*
  printf("\n Eigvals+Eigvects:\n");
  for(int i=0;i<m;i++) 
    {
      printf("%12.9f: (",w[i]);
      for(int j=0;j<n-1;j++) printf("%f, ",z[i*n+j]);
      printf("%f)",z[i*n+n-1]);
      printf("\n");
    }
  */
  delete[] work;delete[] iwork;
  return 0;
#endif
#endif
  return -1;
}



