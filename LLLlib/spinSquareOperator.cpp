#include<spinSquareOperator.hpp>

/*
static members
*/
const double spinSquareOperator::MIN_CONTRIB=1.e-8;   
const double spinSquareOperator::MIN_CONTRIB_SINGLE=1.e-5;   

spinSquareOperator::spinSquareOperator(Basis *new_basis) : Operator(new_basis,0)
{
  doNotTestEqSpat=1;
  ctorBody();
}

spinSquareOperator::spinSquareOperator(Basis *new_basis, 
		    int new_type) : Operator(new_basis,new_type)
{
  doNotTestEqSpat=0;
  ctorBody();  
}

void spinSquareOperator::ctorBody()
{
  /* This operator can handle only (linear combinations of)
     yosStates, i.e. endBasis has to be of type yosBasis. */
  if(endBasis->getBasisType()!=YOS_BASIS) 
  {
  	CxErrors(__FILE__, __LINE__);
  }
  endYosBasis=(yosBasis *)endBasis;
  Nm=endYosBasis->getNm();
  Ne=endYosBasis->getNe();
}

/* Computes the expectation value of S^2 for the i_st-th vector 
   of my_basis. Can deal with non-normalized vectors too. */
double spinSquareOperator::totSpin(int i_st)
{
  double totSsq,norm;
  std::complex<double> a;
  matEl_nonzero(i_st,i_st,&a);
  norm=my_basis->norm(i_st);
  if(norm==0) return 0.0;  /* you tried to compute expectation
			      value of a zero vector; anyway
			      this is strange as basis should
			      not contain zero vectors... */
  totSsq=a.real()/norm/norm;  // 'Normalize' the vector
  return (sqrt(1+4*totSsq)-1)/2.0; // compute S from S(S+1)  
}

/* Computes the expectation value of S^2 for the vector 
   given by coefficients *inp_st with respect to endYosBasis.
   Can deal with non-normalized vectors too. 
   For the sake of security, dimension of *inp_st has to be
   given (which must match the dimension of endYosBasis.
*/
double spinSquareOperator::totSpin(double *inp_st,int check_dim)
{
  double totSsq=0;   /* totalS^2, i.e. S(S+1) */
  long int n=endYosBasis->dimension();
  double *st;
  double prod;
  std::complex<double> a(0,0);

  
  if(n!=check_dim) {
  	CxErrors(__FILE__, __LINE__);
  }
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
	  totSsq = totSsq+prod*a.real();
	}
    }
  delete[] st;
  return (sqrt(1+4*totSsq)-1)/2.0; // compute S from S(S+1)
  //return 0.5;
}


/* Overloaded virual method of Operator... */
int spinSquareOperator::matEl_nonzero(int i_st1,
				      int i_st2,
				      std::complex<double> *a)
{
  double ME=0;   /* the matrix element */

  (*a)=std::complex<double>(0,0);
  if(!spinMatEl(i_st1,i_st2,&ME)) return 0;
  (*a)+=ME;
  return 1;
}


/* The actual calculation of the S^2 matrix element: between i_st1-th and
   i_st2-th state of endYosBasis. The result is put to *ME and the return 
   value is 1 if the MatEl was non-zero and 0 if it was zero.
*/
int spinSquareOperator::spinMatEl(int i_st1,int i_st2,double *ME)
{
  int NrDiffSpins,s1i,s1j,s2i,s2j;

  (*ME)=0;
  st1=(*endYosBasis)[i_st1];st2=(*endYosBasis)[i_st2];  
  if(!doNotTestEqSpat)
    if(!equalSpat(st1,st2)) return 0; // spatial parts of st1,st2 have to be equal
  NrDiffSpins=0;
  for(int k=0;k<Ne;k++)             // find the differences in the spin parts
    if(st1->spin(k)!=st2->spin(k))  
      NrDiffSpins++;
  if(NrDiffSpins==0)
    (*ME)+=3.*Ne/4.;           // Diagonal term.
  if(NrDiffSpins>2) return 0;       // it may not match at max. 2 places
  for(int i=0;i<Ne;i++)
    {	  
      for(int j=0;j<Ne;j++)  // S_i*S_i=3/4 already included, i.e. start with j=i+1
	{
	  if(i==j) continue;
	  s1i=st1->spin(i);s1j=st1->spin(j);s2i=st2->spin(i);s2j=st2->spin(j);
	  /* The spin parts may differ at i-th and j-th position only */
	  if(NrDiffSpins==2)
	    if(!((s1i!=s2i) && (s1j!=s2j))) continue;
	  if(NrDiffSpins==1)
	    if(!((s1i!=s2i) || (s1j!=s2j))) continue;
	  if(st1->j(i)==st1->j(j)) 
	    { // BEWARE! Assume, that only |i+i-...> may exist, not |i+i+...>
	      // which is not allowed due to the Pauli principle
		    if((s1i==s2i) && (s1j==s2j))		  
		      {(*ME)+=-3./4.;continue;}
		    if((s1j==s2i) && (s1i==s2j))		  
		      {(*ME)+=3./4.;continue;}
	    }
	  if((s1i==1) && (s1j==0) && (s2i==0) && (s2j==1))
	    (*ME)+=1.;				
	  if((s1i==1) && (s1j==1) && (s2i==1) && (s2j==1))
	    (*ME)+=1./4.;				
	  if((s1i==0) && (s1j==0) && (s2i==0) && (s2j==0))
	    (*ME)+=1./4.;				
	  if((s1i==0) && (s1j==1) && (s2i==0) && (s2j==1))
	    (*ME)+=-1./2.;
	}
      
    }
  return 1;
}




/* Normalize the n-dimensional state in st_in and put the
   result into st_out. Returns 1 if the state is the zero 
   vector, otherwise returns 0. */
int spinSquareOperator::normalize(double *st_out,double *st_in,int n)
{
  double x=0;
  for(int i=0;i<n;x+=st_in[i]*st_in[i],i++);
  if(x==0) {for(int i=0;i<n;st_out[i]=0,i++);return(1);}
  for(int i=0;i<n;st_out[i]=st_in[i]/sqrt(x),i++);  
  return(0);
}



















