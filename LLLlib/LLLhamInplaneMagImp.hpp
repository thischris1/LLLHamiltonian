/*! \file LLLhamMagImp.hpp */
/* Class for representing LLLhamiltonians with a magnetic impurity of type

   H_imp=\mu*B_x*\sum_{i=1}^{N_e} g(x_i)*(S^+_i+S^-_i) +
         \mu*B*g*\sum_{i=1}^{N_e} S^z_i,

   g(x_i)=g_1*cos(2*\pi*x_i/a),

   i.e. strong (constant) perpendicular field (second, Zeeman term) and small 
   spatially varying magnetic field along the x-direction (S^x=1/2(S^{+}+S^-).

 */

#ifndef LLLhamInplaneMagImp_HPP
#define LLLhamInplaneMagImp_HPP

#include<Basis.hpp>
#include<yosBasis.hpp>
#include<Operator.hpp>
#include<LLLhamiltonian.hpp>

class LLLhamInplaneMagImp : public LLLhamiltonian
{
private:
  double MIstrength;     // 
  double IMIstrength;     // 
  double ZeemanStrength;

public:
  LLLhamInplaneMagImp(yosBasis *new_basis,int new_type,
		      double bli_new,double a_new, double b_new,
		      double new_ZeemanStrength,
		      double new_MIstrength,
		      double new_IMIstrength);

  ~LLLhamInplaneMagImp();

  int set_IMIstrength(double new_IMIstrength);
  double get_IMIstrength();

  int set_MIstrength(double new_MIstrength);
  double get_MIstrength();

  int set_ZeemanStrength (double new_ZeemanStrength);
  double get_ZeemanStrength();

  int matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl);

};

#endif





