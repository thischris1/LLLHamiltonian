#include<LLLhamMagImp.hpp>

LLLhamMagImp::LLLhamMagImp(yosBasis *new_basis,int new_type,double bli_new,
			   double a_new, double b_new,double new_MIstrength) : 
  LLLhamiltonian(new_basis,new_type,bli_new,a_new,b_new)
{
  set_MIstrength(new_MIstrength);
}

LLLhamMagImp::~LLLhamMagImp() {}

int LLLhamMagImp::set_MIstrength(double new_MIstrength)
{
  MIstrength=new_MIstrength;
  return 0;
}

double LLLhamMagImp::get_MIstrength()
{return MIstrength;}

int LLLhamMagImp::matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl)
{
  double MEreal,MEimag,tmp;
  //(*matEl)=std::complex<double>(0,0);
  if(i_st1!=i_st2) return 0;
  MEreal=0;MEimag=0;
  st1=(*LLL_basis)[i_st1];
  for(int i=0;i<Ne;i++)
    {
      tmp=exp(-M_PI/2/Nm/(a/b))*cos(2*M_PI*(st1->j(i)*1.)/(1.*Nm));
      if(st1->spin(i))
	MEreal+=tmp;
      else 
	MEreal-=tmp;
    }
  (*matEl)+=MIstrength*MEreal;
  return 1; // This matrix element is non-zero
}

