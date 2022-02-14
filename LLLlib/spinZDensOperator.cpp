#include <spinZDensOperator.hpp>

spinZDensOperator::spinZDensOperator(const Basis &newBasis,
				     const eigSt &state_to_eval,int new_type) 
  : DensOperator(newBasis,state_to_eval,new_type,true)
{
  spinZOp=new spinZOperator(endYosBasis,new_type);
  setupLandauMatrixForState();
}
 
spinZDensOperator::~spinZDensOperator()
{
  delete spinZOp;
}

std::complex<double> spinZDensOperator::getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2)
{
  std::complex<double> a(0,0);
  spinZOp->matEl_nonzero(i1,i2,&a);
  return std::complex<double>(a.real(), 0.0);
}


int spinZDensOperator::whatTypeOfLandauMatrixDoINeed()
{ return 5;}


