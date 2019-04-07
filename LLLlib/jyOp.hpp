/*!
  \file jyOp.hpp
  Declaration of current density in y-direction
*/

#ifndef JYOP_HPP
#define JYOP_HPP

#include <OneParticleOperator.hpp>

/*!
  \class jyOp : public OneParticleOperator
  \brief class representing current density operator in y-direction
*/
class jyOp : public OneParticleOperator
{

private:
  /*!
    \brief aspect ration
  */
  double m_aTob;

  /*!
    \brief Fluxes of solenoids 1,2 in units of h/e    
   */
  double m_dAlpha1, m_dAlpha2;


public:
  /*!
    \brief constructor
    basis_to_use : basis to be used to calculate expectation value of operator in.
    state_to_eval: state to be evaluated
    aTob: aspect ratio of unit cell
    mat_type: allways 0 (until now)
  */
  jyOp ( const Basis &basis_to_use, 
	 const eigSt &state_to_eval, 
	 double aTob, 
	 int mat_type );

  /*!
    \brief constructor (not yet testet)
    needed to use landau-matrix of ancestor
    basis_to_use : basis to be used to calculate expectation value of operator in.
    ancestor: another OneParticleOperator-derived operator for the state to be evaluated
    aTob: aspect ratio of unit cell
    mat_type: allways 0 (until now)
  */
  jyOp ( const Basis &basis_to_use, 
	 const OneParticleOperator &ancestor, 
	 double aTob,
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
  ~jyOp ();

  /*!
    \brief overides getOneParticleEl of OneParticleOperator to implement
    y-component of current density operator
  */
  std::complex<double> getOneParticleEl (int i, int j);

};

#endif
