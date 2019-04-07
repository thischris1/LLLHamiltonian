/*! \file LLLhamMagImpY.hpp */
/*! \class LLLhamMagImpY

\brief Class for representing LLLhamiltonians with a magnetic impurity of type

   H_imp = \mu*g*\sum_{i=1}^{N_e} B_x(y_i)*(S^+_i+S^-_i) +
           \mu*g*\sum_{i=1}^{N_e} B_z(y_i)*S^z_i +
           \mu*B*g*\sum_{i=1}^{N_e} S^z_i,


   g(y_i)=g_1*cos(2*\pi*y_i/b)
   
   i.e. ('strong') homogeneous perpendicular magnetic field B and weak
   fluctuating perpendicular (B_z) and inplane (B_x) magnetic field. 

   The proportionality constants (in units of energy) are 
   ZeemanStrength (B), MIYstrength (B_z) and IMIYstrength (B_x).
*/

#ifndef LLLhamMagImpY_HPP
#define LLLhamMagImpY_HPP

#include<Basis.hpp>
#include<yosBasis.hpp>
#include<Operator.hpp>
#include<LLLhamiltonian.hpp>

class LLLhamMagImpY : public LLLhamiltonian
{
private:
  double MIYstrength;     // 
  double IMIYstrength;     // 
  double ZeemanStrength;

public:
  LLLhamMagImpY(yosBasis *new_basis,int new_type,
		double bli_new,double a_new, double b_new,
		double new_ZeemanStrength,
		double new_MIYstrength,
		double new_IMIYstrength);
    
  ~LLLhamMagImpY();

  int set_IMIYstrength(double new_IMIYstrength);
  double get_IMIYstrength();

  int set_MIYstrength(double new_MIYstrength);
  double get_MIYstrength();

  int set_ZeemanStrength (double new_ZeemanStrength);
  double get_ZeemanStrength();

  int matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *matEl);

};

#endif

