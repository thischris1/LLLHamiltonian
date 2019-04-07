#include <spinDensOperator.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>
spinDensOperator::spinDensOperator(const Basis &newBasis,
				   const eigSt &state_to_eval,int new_type) 
  : DensOperator(newBasis,state_to_eval,new_type,true)
{
  spSqOp=new spinSquareOperator(endYosBasis);
  setupLandauMatrixForState();
}
 
spinDensOperator::~spinDensOperator()
{
  delete spSqOp;
}

std::complex<double> spinDensOperator::getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2)
{
  std::complex<double> a(0,0);
  //int shuffled=0;
  //if(i1!=i2) shuffled=endYosBasis->shuffleSpins(i1,k1,k2);
  spSqOp->matEl_nonzero(i1,i2,&a);
  //if(shuffled) endYosBasis->restoreSpins(i1);
  return a;
}

int spinDensOperator::whatTypeOfLandauMatrixDoINeed()
{ return 2;}


/* For the purpose of spinDensOperator the states should be 
   compared only by the j part (i.e. without taking spin into
   account. Spins are treated separately (see
   getOneElementOfLandauMatrix(). */

inline int spinDensOperator::s1_nr(int i) const
{return (int) st1->j(i)*2;}
inline int spinDensOperator::s2_nr(int i) const
{return (int) st2->j(i)*2;} 
