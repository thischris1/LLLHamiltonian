#include<LLLhamMagImpY.hpp>

LLLhamMagImpY::LLLhamMagImpY(yosBasis *new_basis,int new_type,
			     double bli_new,double a_new, double b_new,
			     double new_ZeemanStrength,
			     double new_MIYstrength,
			     double new_IMIYstrength) : 
  LLLhamiltonian(new_basis,new_type,bli_new,a_new,b_new)
{
  set_MIYstrength(new_MIYstrength);
  set_IMIYstrength(new_IMIYstrength);
  set_ZeemanStrength(new_ZeemanStrength);
}

LLLhamMagImpY::~LLLhamMagImpY() {}

int LLLhamMagImpY::set_IMIYstrength(double new_IMIYstrength)
{
  IMIYstrength=new_IMIYstrength;
  return 0;
}
double LLLhamMagImpY::get_IMIYstrength()
{return IMIYstrength;}

int LLLhamMagImpY::set_MIYstrength(double new_MIYstrength)
{
  MIYstrength=new_MIYstrength;
  return 0;
}
double LLLhamMagImpY::get_MIYstrength()
{return MIYstrength;}

int LLLhamMagImpY::set_ZeemanStrength (double new_ZeemanStrength)
{
  ZeemanStrength=new_ZeemanStrength;
  return 0;
}

double LLLhamMagImpY::get_ZeemanStrength()
{return ZeemanStrength;}

int LLLhamMagImpY::matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl)
{
  double MEreal,MEimag,tmp;
  int nr_diff,idiff=0,itmp;   // for the B_x inhmgty

  //(*matEl)=std::complex<double>(0,0);
  MEreal=0;MEimag=0;
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
      return 1;
    }
  nr_diff=0; // in how many spins do st1,st2 differ
  for(int i=0;i<Ne;i++)  // the inhomogeneous B_x term
    if(st1->s1_nr(i)!=st2->s1_nr(i)) {idiff=i;nr_diff++;}
  if(nr_diff>1) return 0;  // inhmgty is a 1-particle operator
  if(nr_diff==0) return 0; // actually, this should not happen (i_st1!=i_st2)
  itmp=st1->j(idiff)-st2->j(idiff);
  if(itmp*itmp!=1) return 0; // they should differ in one j: j'=j+1 or j-1
  itmp=st1->spin(idiff)+st2->spin(idiff);
  if(itmp==0) tmp=IMIYstrength;  // spins at j,j' are different
  else
      tmp=MIYstrength*itmp/2;    // the spins are equal
  (*matEl)+=tmp;
  return 1; // This matrix element is non-zero
}
