/*!
  \file spinXDensOperator.hpp
*/
#ifndef SPINXDENSOPERATOR_HPP
#define SPINXDENSOPERATOR_HPP

#include <densOperator.hpp>
#include <spinXOperator.hpp>
#include <fstream>
#include <iostream>

class spinXDensOperator : public DensOperator
{
private:
  spinXOperator *spinXOp;
protected:
  virtual std::complex<double> getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2);
  virtual int whatTypeOfLandauMatrixDoINeed();
public:
  spinXDensOperator(const Basis &newBasis,const eigSt &state_to_eval,int new_type);
  virtual ~spinXDensOperator();
};

#endif
