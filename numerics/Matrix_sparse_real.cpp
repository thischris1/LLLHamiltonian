/*!
  \file Matrix_sparse_real.cpp
 \brief Implementation of Class Matrix_sparse_real for sparse matrix diagonalization
 \author Karel Vyborny
*/
#include <numerics/Matrix_sparse_real.hpp>
#include <stdio.h>
#include <LLLlib/LLLlib.h>


/*
  static members
 */
const double Matrix_sparse_real::DIAG_TOLERANCE=1.e-8; // precision at diagonalization


/*************************************
 Functions related to the class (and to f02fjf_ of the NAG library)
**************************************/

// !!!!!!!!! See comments in the .hpp file !!!!!!!!!!!!

void GoToImage(int *iflag,int *n,double z[],double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
{
  ((Matrix_sparse_real *) iwork)->image(iflag,z,w);
}


// Compute scalar product of two vectors 

double GoToDot(int *iflag,int *n,double z[],double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
{
  double sum=0;
  for(int i=0;i<*n;i++) sum+=z[i]*w[i];
  return sum;
  // A slower alternative...
  //return ((Matrix_sparse_real *) iwork)->dot(iflag,z,w);
}


void GoToMonit(int *istate,int *nextit,int *nevals,int *nevecs,int *k,double f[],double d[])
{
  return;

}


/**********************************************
 Methods
***********************************************/

Matrix_sparse_real::Matrix_sparse_real(int *new_row_NrNz,int *new_col_ind,double *new_mat_sp,int new_nz_elements,int new_n,int new_how_many_vects,int add_to_workspace)
{
  row_NrNz=new_row_NrNz;
  col_ind=new_col_ind;
  mat_sp=new_mat_sp;
  nz_elements=new_nz_elements;
  dimension=new_n;
  how_many_eigvects=new_how_many_vects;
  add_to_diag_workspace_sparse=add_to_workspace;
}

Matrix_sparse_real::~Matrix_sparse_real()
{}

void Matrix_sparse_real::set_k(int k)
{
  add_to_diag_workspace_sparse=k;
}

int Matrix_sparse_real::foomax(int a,int b)
{
  if(a>b) return a;
  else return b;
}


/*! 
Compute y=A*x; a slight trouble is that only one half (upper, say) of the matrix
   is stored (say A=U+D+L=U+D+U^T; stored is just D,U). Thus, to get Ax you have
   to compute U*x and U^T*x.
   Written according to
http://www.netlib.org/linalg/html_templates/node98.html#SECTION00932100000000000000
*/
void Matrix_sparse_real::image(int *iflag,double z[],double w[])
{
  int i_Nz=0;
  double matEl;
  int cInd,elsInRow;

  for(int i=0;i<dimension;i++) w[i]=0;  

  for(int i=0;i<dimension;i++)
    {
      //w[i]=0;
      elsInRow=row_NrNz[i];
      /* D*x would be counted twice, once here and once at U^T*x, otherwise */    
      if(col_ind[i_Nz]==i) w[i]-=z[i]*mat_sp[i_Nz]; 
      for(int j=i_Nz;j<i_Nz+elsInRow;j++)
	{
	  matEl=mat_sp[j];
	  cInd=col_ind[j];
	  w[i]+=z[cInd]*matEl; // U*x
	  w[cInd]+=z[i]*matEl; // U^T*x
	}
      i_Nz+=elsInRow;
    }
}


/* Compute scalar product of two vectors 
   Not used at the moment: GoToDot() is enough as it does not need to 
   know anything about the matrix.
*/
double Matrix_sparse_real::dot(int *iflag,double z[],double w[])
{
  double sum=0;
  for(int i=0;i<dimension;i++) sum+=z[i]*w[i];
  return sum;
}

/* Can be used to monitor progress of the diagonalization
   Also not used at the moment...
*/
void Matrix_sparse_real::monit(int *istate,int *nextit,int *nevals,int *nevecs,int *k,double f[],double d[])
{
  return;
}


/* Guess what... parameters: 
   x       eigenvectors, stored one after another in the 1D array
           (size has to be at least (k+4)*dimension???).
   eigval  eigenvalues
*/
int Matrix_sparse_real::diagonalize(double *x,int novecs,double *eigval)
{
	#ifdef NAG
	glLogger.debug("Entering Matrix_sparse_real::diagonalize");
  int n,m,k,noits,nrx,lwork,lrwork,liwork,ifail;
  int *iwork;
  double tol;
  double *work,*rwork;
  
// Different implementation for Intel MKL
  // See documentation of the NAG routine f02fjf for explanation of the
  // variables...

  n=dimension;         
  m=how_many_eigvects;
  k=m+add_to_diag_workspace_sparse;
  if(k>n) k=n;
  noits=DIAG_MAX_ITERATIONS;
  tol=DIAG_TOLERANCE;
  //novecs=0;
  nrx=n;
  //x=new double[nrx*k];
  lwork=3*k+foomax(k*k,2*n);
  work=allocK<double>(lwork);
  lrwork=1;
  rwork=new double[lrwork];
  liwork=1;
  //iwork=new int[liwork];
  iwork=(int *)this;
  ifail=-1;                  // Soft noisy exit from the diag. if there's error
  //printf("n:%d, m:%d, k:%d, nrx:%d, lwork:%d, lrwork:%d, liwork:%d.\n",n,m,k,nrx,lwork,lrwork,liwork);

  f02fjf_(&n,&m,&k,&noits,&tol,GoToDot,GoToImage,GoToMonit,&novecs,x,&nrx,eigval,work,&lwork,rwork,&lrwork,iwork,&liwork,&ifail);

  glLogger.info(", exiting from f02fjf()");
 
  delete[] work;
  return ifail;
#endif
  throw (CxErrors("Not implemented for this configuration"));
  return -1;
}







