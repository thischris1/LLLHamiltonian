#include <jyOp.hpp>

jyOp::jyOp ( const Basis &basis_to_use, 
	     const eigSt &state_to_eval,
	     double aTob,
	     int mat_type ) :
  OneParticleOperator (basis_to_use, state_to_eval, mat_type, true, false)
{
  m_aTob = aTob;
}
 
jyOp::jyOp ( const Basis &basis_to_use, 
	     const OneParticleOperator &ancestor,
	     double aTob,
	     int mat_type ) :
  OneParticleOperator (basis_to_use, ancestor, mat_type, true)
{
  m_aTob = aTob;
}

void jyOp::setSolenoidFluxes (double alpha1, double alpha2)
{
  m_dAlpha1 = alpha1;
  m_dAlpha2 = alpha2;
}


jyOp::~jyOp ()
{
}

std::complex<double> jyOp::getOneParticleEl (int i, int j)
{
  return std::complex<double>(0.0, 0.5) * 
    ( endYosBasis->getCachedWF(j) * conj( endYosBasis->getCachedWFdy(i) ) 
      - conj(  endYosBasis->getCachedWF(i) ) * endYosBasis->getCachedWFdy(j) )
    + m_aTob * (m_dAlpha2 - 2.0*M_PI*Nm*x) * conj( endYosBasis->getCachedWF(i) ) * endYosBasis->getCachedWF(j);
}
