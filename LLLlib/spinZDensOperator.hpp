/*!
  \file spinZDensOperator.hpp
*/
#ifndef SPINZDENSOPERATOR_HPP
#define SPINZDENSOPERATOR_HPP

#include <densOperator.hpp>
#include <spinZOperator.hpp>
#include <fstream>
#include <iostream>

class spinZDensOperator : public DensOperator
{
private:
  spinZOperator *spinZOp;
protected:
  virtual std::complex<double> getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2);
  virtual int whatTypeOfLandauMatrixDoINeed();
public:
  spinZDensOperator(const Basis &newBasis,const eigSt &state_to_eval,int new_type);
  virtual ~spinZDensOperator();
};

#endif
