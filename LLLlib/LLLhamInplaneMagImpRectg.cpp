#include<LLLhamInplaneMagImpRectg.hpp>

LLLhamInplaneMagImpRectg::LLLhamInplaneMagImpRectg(
				   yosBasis *new_basis,int new_type,
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

LLLhamInplaneMagImpRectg::~LLLhamInplaneMagImpRectg() {}

int LLLhamInplaneMagImpRectg::set_IMIstrength(double new_IMIstrength)
{
  IMIstrength=new_IMIstrength;
  return 0;
}
double LLLhamInplaneMagImpRectg::get_IMIstrength()
{return IMIstrength;}

int LLLhamInplaneMagImpRectg::set_MIstrength(double new_MIstrength)
{
  MIstrength=new_MIstrength;
  return 0;
}
double LLLhamInplaneMagImpRectg::get_MIstrength()
{return MIstrength;}

int LLLhamInplaneMagImpRectg::set_ZeemanStrength (double new_ZeemanStrength)
{
  ZeemanStrength=new_ZeemanStrength;
  return 0;
}
double LLLhamInplaneMagImpRectg::get_ZeemanStrength()
{return ZeemanStrength;}

int LLLhamInplaneMagImpRectg::matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl)
{
  double MEreal,MEimag,tmp;
  int nr_diff,jdiff=0;   // for the B_x inhmgty
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
      for(int i=0;i<Ne;i++)  // the inhomogeneous B_z term (also only diagonal)
	{
	  tmp=exp(-M_PI/2/Nm/(a/b));
	  if(inhmgtyGeometry(st1->j(i)))
	    tmp=-tmp;
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
  tmp=exp(-M_PI/2/Nm/(a/b));
  //if(inhmgtyGeometry((jdiff) % Nm))
  if(inhmgtyGeometry((jdiff+1) % Nm))
    tmp=-tmp;
  MEreal+=tmp;
  (*matEl)+=IMIstrength*MEreal;
  return 1; // This matrix element is non-zero
}

inline int LLLhamInplaneMagImpRectg::inhmgtyGeometry(int j)
{
  /* rectg
     non-centered, B=+ in the middle, middle bigger;
     i.e. (at Ne/Nm=6/9): B=- at j=0,1,7,8, B=+ at j=2,3,4,5,6 */
  //if(j<Nm/4 || j>3*(Nm/4)) return 1;

  /* rectg2
     centered, B=+ at the edges, middle smaller;
     i.e. (at Ne/Nm=6/9): B=- at j=3,4,5,6 B=+ at j=0,1,2,7,8, or
     (at 8/12): B=- at j=4,5,6,7,8,9, B=+ at j=0,1,2,3,10,11 */
  if(j>Nm/4 && j<=3*(Nm/4)) return 1;

  /* rectg3
     centered, B=+ in the middle, edges smaller;
     i.e. (at Ne/Nm=6/9): B=+ at j=3,4,5,6 B=- at j=0,1,2,7,8 */
  //if(!(j>Nm/4 && j<=3*(Nm/4))) return 1;

  /* rectg4
     centered, B=+ at the edges, middle smaller; domain sizes 2:1
     i.e. (at Ne/Nm=6/9): B=- at j=4,5,6 B=+ at j=0,1,2,3,7,8
     (at 8/12): B=- at j=5,6,7,8, B=+ at j=0,1,2,3,4,9,10,11 */
  //if(j>Nm/3 && j<=2*(Nm/3)) return 1;

  return 0;
}





