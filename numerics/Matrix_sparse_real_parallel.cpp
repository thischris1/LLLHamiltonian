// Class for sparse matrix diagonalization

#ifdef COMPILE_PARALLELIZED

#include<Matrix_sparse_real_parallel.hpp>
#include<stdio.h>
#include<myaccessories.hpp>


/*
  static members
 */
const double Matrix_sparse_real_parallel::DIAG_TOLERANCE=1.e-8; 
                                            // precision at diagonalization


/*************************************
 Functions related to the class (and to f02fjf_ of the NAG library)
**************************************/

// !!!!!!!!! See comments in the .hpp file !!!!!!!!!!!!

void GoToImage_parallel(int *iflag,int *n,double z[],double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
{
  ((Matrix_sparse_real_parallel *) iwork)->image(iflag,z,w);
}


// Compute scalar product of two vectors 

double GoToDot_parallel(int *iflag,int *n,double z[],double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
{
  double sum=0;
  for(int i=0;i<*n;i++) sum+=z[i]*w[i];
  return sum;
  // A slower alternative...
  //return ((Matrix_sparse_real *) iwork)->dot(iflag,z,w);
}

// Can be used to monitor progress of the diagonalization

void GoToMonit_parallel(int *istate,int *nextit,int *nevals,int *nevecs,int *k,double f[],double d[])
{
  return;
  // A slower alternative...
  //((Matrix_sparse_real *) iwork)->monit(istate,nextit,nevals,nevecs,k,f,d);
  //return 0;
  //printf("Step %d finished. Eigvals so far: %d, eigvecs:%d.\n",*nextit,*nevals,*nevecs);
  //fflush(stdout);
  //return NULL;
}



/**********************************************
 Methods
***********************************************/

Matrix_sparse_real_parallel::Matrix_sparse_real_parallel(int *new_p_row_NrNz,
			      std::vector<MatEl> **new_p_rowOfMatEls,
			      int new_nz_elements,
			      int new_dimension,
			      int new_rowEvaluated,
			      int new_eigsToFind,
			      int new_add_to_diag_workspace_sparse)
{
  p_row_NrNz=new_p_row_NrNz;
  p_rowOfMatEls=new_p_rowOfMatEls;
  nz_elements=new_nz_elements;
  dimension=new_dimension;
  rowEvaluated=new_rowEvaluated;
  eigsToFind=new_eigsToFind;
  add_to_diag_workspace_sparse=new_add_to_diag_workspace_sparse;

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  //req_productReady = new MPI_Request[numprocs];
  if(myid==0) vecTmp = allocK<double>(dimension);
}

Matrix_sparse_real_parallel::~Matrix_sparse_real_parallel()
{
  //delete []req_productReady;
  if(myid==0) delete []vecTmp;
}

void Matrix_sparse_real_parallel::set_k(int k)
{
  add_to_diag_workspace_sparse=k;
}

int Matrix_sparse_real_parallel::foomax(int a,int b)
{
  if(a>b) return a;
  else return b;
}


/* Compute y=A*x; a slight trouble is that only one half (upper, say) of the matrix
   is stored (say A=U+D+L=U+D+U^T; stored is just D,U). Thus, to get Ax you have
   to compute U*x and U^T*x.
   Written according to
http://www.netlib.org/linalg/html_templates/node98.html#SECTION00932100000000000000
*/
void Matrix_sparse_real_parallel::image(int *iflag,double z[],double w[])
{

  // Summation within the master proc
  image_part(z,w);

  // Collect results from slave procs

  //int allProductsReady = 0;
  MPI_Status tmp_status;
  for(int i=1;i<numprocs;i++)
    {
      MPI_Recv(vecTmp,dimension,MPI_DOUBLE,i,
		MPICODE_PRODUCT_READY,MPI_COMM_WORLD,&tmp_status);
      for(int j=0;j<dimension;j++)
	z[j]+=vecTmp[j]; // ... and sum them up
    }
}


void Matrix_sparse_real_parallel::image_part(double vecIn[],double vecOut[])
{
  double matEl;
  int cInd,elsInRow;
  std::vector<MatEl> *p_thisRow;

  for(int i=0;i<dimension;i++) vecOut[i]=0;  

  int currRow = myid; // Nr. of the row in reality
  for(int i=0;i<rowEvaluated;i++)
    {
      elsInRow = p_row_NrNz[i];
      p_thisRow = p_rowOfMatEls[i];
      /* D*x would be counted twice, once here and once at U^T*x, otherwise */    
      if(((*p_thisRow)[0]).col==currRow) 
	vecOut[currRow] -=vecIn[currRow]*((*p_thisRow)[0]).val; 
      for(int j=0;j<elsInRow;j++)
	{
	  matEl = ((*p_thisRow)[j]).val;
	  cInd  = ((*p_thisRow)[j]).col;
	  vecOut[currRow] += vecIn[cInd]*matEl; // U*x
	  vecOut[cInd] += vecIn[currRow]*matEl; // U^T*x
	}
      currRow += numprocs;
      if(currRow >= dimension) 
	{
	  glLogger.error("Matrix_sparse_real_parallel::image_part");
	  glLogger.error("currRow > dimension.");
	  exit(1);
	}
    }
}


/* Compute scalar product of two vectors 
   Not used at the moment: GoToDot() is enough as it does not need to 
   know anything about the matrix.
*/
double Matrix_sparse_real_parallel::dot(int *iflag,double z[],double w[])
{
  double sum=0;
  for(int i=0;i<dimension;i++) sum+=z[i]*w[i];
  return sum;
}

/* Can be used to monitor progress of the diagonalization
   Also not used at the moment...
*/
void Matrix_sparse_real_parallel::monit(int *istate,int *nextit,int *nevals,
					int *nevecs,int *k,double f[],double d[])
{
  return;
}


/* Guess what... parameters: 
   x       eigenvectors, stored one after another in the 1D array
           (size has to be at least (k+4)*dimension???).
   eigval  eigenvalues
*/
int Matrix_sparse_real_parallel::diagonalize(double *x,int novecs,double *eigval)
{
  if(myid>0) // Slave process
    {
      int whatToDo;
      double *vecIn=allocK<double>(dimension);
      double *vecOut=allocK<double>(dimension);
      MPI_Status my_status;
      // Wait for request to compute A*vecIn=vecOut
      MPI_Recv(&whatToDo,1,MPI_INT,0,
	       MPICODE_SLAVE_WHATTODO,MPI_COMM_WORLD,&my_status);
      while(whatToDo>0)  // whatToDo=0 ... diagonalization finished
	{
	  // Get vecIn
	  MPI_Recv(vecIn,dimension,MPI_DOUBLE,0,
		   MPICODE_COMPUTE_PRODUCT,MPI_COMM_WORLD,&my_status);
	  // Compute A*vecIn
	  image_part(vecIn,vecOut);
	  // Send the result back to proc 0
	  MPI_Ssend(vecOut,dimension,MPI_DOUBLE,0,
		    MPICODE_PRODUCT_READY,MPI_COMM_WORLD);
	  // and wait for next orders
	  MPI_Recv(&whatToDo,1,MPI_INT,0,
		   MPICODE_SLAVE_WHATTODO,MPI_COMM_WORLD,&my_status);
	}
      return 0; // Slave procs finish here
    }
  
  // Master process

  int n,m,k,noits,nrx,lwork,lrwork,liwork,ifail;
  int *iwork;
  double tol;
  double *work,*rwork;
 
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
  //printf("n:%d, m:%d, k:%d, nrx:%d, lwork:%d, lrwork:%d, liwork:%d.\n",
  //	 n,m,k,nrx,lwork,lrwork,liwork);
  f02fjf_(&n,&m,&k,&noits,&tol,
	  GoToDot_parallel,GoToImage_parallel,GoToMonit_parallel,
	  &novecs,x,&nrx,eigval,work,&lwork,rwork,&lrwork,iwork,&liwork,&ifail);
  glLogger.info(", exiting from f02fjf()");
  /*
  printf("\n Exit code:%d.\n",ifail);
  printf("\n Eigvals+Eigvects:\n");
  for(int i=0;i<m;i++) 
    {
      printf("%f: (",eigval[i]);
      for(int j=0;j<n-1;j++) printf("%f, ",x[i*n+j]);
      printf("%f)",x[i*n+n-1]);
      printf("\n");
    } 
  */
  delete[] work;
  return ifail;
}


#endif




