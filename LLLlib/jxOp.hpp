/*!
  \file OneParticleOperator.hpp
*/

#ifndef JXOP_HPP
#define JXOP_HPP

#include <OneParticleOperator.hpp>

/*!
  \class jxOp : public OneParticleOperator
  \brief class representing current density operator in x-direction
*/
class jxOp : public OneParticleOperator
{

private:

  /*!
    \brief Fluxes of solenoids 1,2 in units of h/e    
   */
  double m_dAlpha1, m_dAlpha2;


public:
  /*!
   \fn jxOp ( const Basis &basis_to_use, const eigSt &state_to_eval, int mat_type )
   \brief constructor
   \param basis_to_use : basis to be used to calculate expectation value of operator in.
   \param state_to_eval: state to be evaluated
   \param mat_type: allways 0 (until now)
  */
  jxOp ( const Basis &basis_to_use, 
	 const eigSt &state_to_eval, 
	 int mat_type );


  /*!
  \fn jxOp ( const Basis &basis_to_use, const OneParticleOperator &ancestor, int mat_type)
  \brief constructor (not yet testet) needed to use landau-matrix of ancestor
  \param basis_to_use : basis to be used to calculate expectation value of operator in.
  \param ancestor: another OneParticleOperator-derived operator for the state to be evaluated
  \param mat_type: allways 0 (until now)
  */
  jxOp ( const Basis &basis_to_use, 
	 const OneParticleOperator &ancestor, 
	 int mat_type);


  /*!
   \fn setSolenoidFluxes (double alpha1, double alpha2)
   \brief Set the fluxes in units of h/e of solenoids modifiing boundary conditions.
   alpha1 affects the x-momentum, alpha2 the y-momentum.
  */
  void setSolenoidFluxes (double alpha1, double alpha2);

  
  /*!
    \brief destructor
  */
  ~jxOp ();

  /*!
    \brief overides getOneParticleEl of OneParticleOperator to implement
    x-component of current density operator
  */
  std::complex<double> getOneParticleEl (int i, int j);

};

#endif
