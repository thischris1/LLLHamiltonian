#include <spinDownDensOperator.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>
spinDownDensOperator::spinDownDensOperator(const Basis &newBasis,
					   const eigSt &state_to_eval,int new_type) 
  : DensOperator(newBasis,state_to_eval,new_type,true)
{  
  setupLandauMatrixForState();
}
 
spinDownDensOperator::~spinDownDensOperator()
{}

std::complex<double> spinDownDensOperator::getOneElementOfLandauMatrix(int k1,int k2,
							 int i1,int i2)
{
  if((s1_nr(k1) % 2) != (s2_nr(k2) % 2)) return 0.0;  // spins equal?
  if((s1_nr(k1) % 2) != 0) return 0.0;                // spin down?
  return 1.0; // yes
}

int spinDownDensOperator::whatTypeOfLandauMatrixDoINeed()
{ return 3;}
