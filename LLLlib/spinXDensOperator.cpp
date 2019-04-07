#include <spinXDensOperator.hpp>

spinXDensOperator::spinXDensOperator(const Basis &newBasis,
				     const eigSt &state_to_eval,int new_type) 
  : DensOperator(newBasis,state_to_eval,new_type,true)
{
  spinXOp=new spinXOperator(endYosBasis);
  setupLandauMatrixForState();
}
 
spinXDensOperator::~spinXDensOperator()
{
  delete spinXOp;
}

std::complex<double> spinXDensOperator::getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2)
{
  std::complex<double> a(0,0);
  spinXOp->matEl_nonzero(i1,i2,&a);
  return a;
}


int spinXDensOperator::whatTypeOfLandauMatrixDoINeed()
{ return 4;}

