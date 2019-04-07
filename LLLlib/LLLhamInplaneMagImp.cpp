#include<LLLhamInplaneMagImp.hpp>

LLLhamInplaneMagImp::LLLhamInplaneMagImp(yosBasis *new_basis,int new_type,
					 double bli_new,double a_new, double b_new,
					 double new_ZeemanStrength,
					 double new_MIstrength,
					 double new_IMIstrength) : 
  LLLhamiltonian(new_basis,new_type,bli_new,a_new,b_new)
{
  set_MIstrength(new_MIstrength);
  set_IMIstrength(new_IMIstrength);
  set_ZeemanStrength(new_ZeemanStrength);
}

LLLhamInplaneMagImp::~LLLhamInplaneMagImp() {}

int LLLhamInplaneMagImp::set_IMIstrength(double new_IMIstrength)
{
  IMIstrength=new_IMIstrength;
  return 0;
}
double LLLhamInplaneMagImp::get_IMIstrength()
{return IMIstrength;}

int LLLhamInplaneMagImp::set_MIstrength(double new_MIstrength)
{
  MIstrength=new_MIstrength;
  return 0;
}
double LLLhamInplaneMagImp::get_MIstrength()
{return MIstrength;}

int LLLhamInplaneMagImp::set_ZeemanStrength (double new_ZeemanStrength)
{
  ZeemanStrength=new_ZeemanStrength;
  return 0;
}

double LLLhamInplaneMagImp::get_ZeemanStrength()
{return ZeemanStrength;}

int LLLhamInplaneMagImp::matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl)
{
  double
      MEreal,
      MEimag,
      tmp;

  int
      nr_diff,
      jdiff=0;   // for the B_x inhmgty

  //(*matEl)=std::complex<double>(0,0);
  MEreal=0;
  MEimag=0;
  st1=(*LLL_basis)[i_st1];
  st2=(*LLL_basis)[i_st2];
  if(i_st1==i_st2) 
    {
      for(int i=0;i<Ne;i++)  // the Zeeman term (only diagonal)
	if(st1->spin(i))
	  MEreal+=1.;
	else
	  MEreal-=1.;
      (*matEl)+=ZeemanStrength*MEreal;
      for(int i=0;i<Ne;i++)  // the inhomogeneous B_z term (also only diagonal)
	{
	  tmp=exp(-M_PI/2/Nm/(a/b))*cos(2*M_PI*(st1->j(i)*1.)/(1.*Nm));
	  if(st1->spin(i))
	    MEreal+=tmp;
	  else 
	    MEreal-=tmp;
	}
      (*matEl)+=MIstrength*MEreal;
      return 1;
    }
  if(!equalSpat(st1,st2)) return 0;
  nr_diff=0; // in how many spins do st1,st2 differ
  for(int i=0;i<Ne;i++)  // the inhomogeneous B_x term
    if(st1->spin(i)!=st2->spin(i)) {jdiff=st1->j(i);nr_diff++;}
  if(nr_diff>1) return 0;  // inhmgty is a 1-particle operator
  if(nr_diff==0) return 0; // actually, this should not happen (i_st1!=i_st2)
  /* i.e. the two states differ only in one spin (of a 1-part. state) and 
     this 1-part. state has j=jdiff */
  tmp=exp(-M_PI/2/Nm/(a/b))*sin(2*M_PI*((jdiff)*1.)/(1.*Nm));
  MEreal+=tmp;
  (*matEl)+=IMIstrength*MEreal;
  return 1; // This matrix element is non-zero
}








