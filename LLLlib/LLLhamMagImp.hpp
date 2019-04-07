/*! \file LLLhamMagImp.hpp */
/* Class for representing LLLhamiltonians with a magnetic impurity of type

   H_imp=\mu*B*\sum_{i=1}^{N_e} g(x_i)S^z_i,

   g(x_i)=g_1*cos(2*\pi*x_i/a)
   
   i.e. spatially varying g-factor (or more or less slightly varying B).
   The proportionality constant (i.e. \mu*B*g_1) is called MIstrength.

 */

#ifndef LLLhamMagImp_HPP
#define LLLhamMagImp_HPP

#include<Basis.hpp>
#include<yosBasis.hpp>
#include<Operator.hpp>
#include<LLLhamiltonian.hpp>

class LLLhamMagImp : public LLLhamiltonian
{
private:
  double MIstrength;   /* strength of the magnetic impurity;
			  MIstrength=\mu_B*B*g_1 */

public:
  LLLhamMagImp(yosBasis *new_basis,int new_type,double bli_new,
	       double a_new, double b_new,double new_MIstrength);
  ~LLLhamMagImp();
  int set_MIstrength(double new_MIstrength);
  double get_MIstrength();
  virtual int matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl);
  
};

#endif
