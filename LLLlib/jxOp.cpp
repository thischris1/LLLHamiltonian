/*!
  \file jxOp.cpp
 */

#include <jxOp.hpp>

jxOp::jxOp ( const Basis &basis_to_use, 
	     const eigSt &state_to_eval, 
	     int mat_type ) :
  OneParticleOperator (basis_to_use, state_to_eval, mat_type, true, false)
{
  m_dAlpha1 = 0.0;
  m_dAlpha2 = 0.0;
}
 
jxOp::jxOp ( const Basis &basis_to_use, 
		const OneParticleOperator &ancestor, 
		int mat_type ) :
  OneParticleOperator (basis_to_use, ancestor, mat_type, true)
{
}


jxOp::~jxOp ()
{
}

void jxOp::setSolenoidFluxes (double alpha1, double alpha2)
{
  m_dAlpha1 = alpha1;
  m_dAlpha2 = alpha2;
}

std::complex<double> jxOp::getOneParticleEl (int i, int j)
{
  return std::complex<double>(0.0, 0.5) * 
    ( endYosBasis->getCachedWF(j) * conj( endYosBasis->getCachedWFdx(i) )
     - conj( endYosBasis->getCachedWF(i) ) * endYosBasis->getCachedWFdx(j) )
    + m_dAlpha1 * conj( endYosBasis->getCachedWF(i) ) * endYosBasis->getCachedWF(j);

}
