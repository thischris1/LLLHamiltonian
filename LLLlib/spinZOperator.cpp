#include<spinZOperator.hpp>

/*
  static members
*/
const double spinZOperator::MIN_CONTRIB=1.e-8;   
const double spinZOperator::MIN_CONTRIB_SINGLE=5.e-4;   

spinZOperator::spinZOperator(Basis *new_basis, 
		    int new_type) : Operator(new_basis,new_type)
{
  doNotTestEqSpat=0;
  ctorBody();  
}

/* // This is not what you need for spinZ operator in yosBasis (unlike
   // spinSquareOperator).
   
spinZOperator::spinZOperator(Basis *new_basis) : Operator(new_basis,0)
{
  std::cerr << "spinZOperator::spinZOperator Attempt to use the class with doNotTestEqSpat=1. Implement this in spinZOperator::totSpinZ first.\n ";exit(1);
  doNotTestEqSpat=1;
  ctorBody();
}
*/

void spinZOperator::ctorBody()
{
  /* This operator can handle only (linear combinations of)
     yosStates, i.e. endBasis has to be of type yosBasis. */
  if(endBasis->getBasisType()!=YOS_BASIS) {throw CxErrors(__FILE__, __LINE__);}
  endYosBasis=(yosBasis *)endBasis;
  Nm=endYosBasis->getNm();
  Ne=endYosBasis->getNe();
}

/* Computes the expectation value of S^2 for the i_st-th vector 
   of my_basis. Can deal with non-normalized vectors too. */
double spinZOperator::totSpinZ(int i_st)
{
  double totSz,norm;
  std::complex<double> a;
  matEl_nonzero(i_st,i_st,&a);
  norm=my_basis->norm(i_st);
  if(norm==0) return 0.0;  /* you tried to compute expectation
			      value of a zero vector; anyway
			      this is strange as basis should
			      not contain zero vectors... */
  totSz=a.real()/norm/norm;  // 'Normalize' the vector
  return totSz; 
}

/* Computes the expectation value of S^2 for the vector 
   given by coefficients *inp_st with respect to endYosBasis.
   Can deal with non-normalized vectors too. 
   For the sake of security, dimension of *inp_st has to be
   given (which must match the dimension of endYosBasis.
*/
double spinZOperator::totSpinZ(double *inp_st,int check_dim)
{
  double totSz=0;   /* totalS^2, i.e. S(S+1) */
  long int n=endYosBasis->dimension();
  double *st;
  double prod;
  std::complex<double> a(0,0);

  
  if(n!=check_dim) {throw CxErrors(__FILE__, __LINE__);}
  st=allocK<double>(n);
  normalize(st,inp_st,n);
  for(int i1=0;i1<n;i1++)    // Each state in the linear combination
    {
      if(fabs(st[i1])<MIN_CONTRIB_SINGLE) continue;
      for(int i2=0;i2<n;i2++)  // left with each state in the LC right...
	{
	  prod=st[i1]*st[i2];
	  if(fabs(prod)<MIN_CONTRIB) continue;
	  if(!matEl_nonzero(i1,i2,&a)) continue;
	  totSz = totSz+prod*a.real();
	}
    }
  delete[] st;
  return totSz; 
}


/* Overloaded virtual method of Operator... */
int spinZOperator::matEl_nonzero(int i_st1,
				 int i_st2,
				 std::complex<double> *a)
{
  double ME=0;   /* the matrix element */

  (*a)=std::complex<double>(0,0);
  if(!spinZMatEl(i_st1,i_st2,&ME)) return 0;
  (*a)+=ME;
  return 1;
}


/* The actual calculation of the S^2 matrix element: between i_st1-th and
   i_st2-th state of endYosBasis. The result is put to *ME and the return 
   value is 1 if the MatEl was non-zero and 0 if it was zero.
*/
int spinZOperator::spinZMatEl(int i_st1,int i_st2,double *ME)
{

  (*ME)=0;
  st1=(*endYosBasis)[i_st1];st2=(*endYosBasis)[i_st2];
  /*
  if(!doNotTestEqSpat)
    if(!equalSpat(st1,st2)) return 0; // spatial parts of st1,st2 have to be equal
  */
  if(i_st1!=i_st2) return 0;
  for(int i=0;i<Ne;i++)
    if(st1->spin(i))
      (*ME)+=0.5;
    else
      (*ME)-=0.5;
  return 1;
}




/* Normalize the n-dimensional state in st_in and put the
   result into st_out. Returns 1 if the state is the zero 
   vector, otherwise returns 0. */
int spinZOperator::normalize(double *st_out,double *st_in,int n)
{
  double x=0;
  for(int i=0;i<n;x+=st_in[i]*st_in[i],i++);
  if(x==0) {for(int i=0;i<n;st_out[i]=0,i++);return(1);}
  for(int i=0;i<n;st_out[i]=st_in[i]/sqrt(x),i++);  
  return(0);
}







