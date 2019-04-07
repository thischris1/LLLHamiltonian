/*!
  \file spinDownDensOperator.hpp
*/
#ifndef SPINDOWNDENSOPERATOR_HPP
#define SPINDOWNDENSOPERATOR_HPP

#include <densOperator.hpp>
#include <fstream>
#include <iostream>

class spinDownDensOperator : public DensOperator
//class spinDownDensOperator : public DensOp
{
protected:
  virtual std::complex<double> getOneElementOfLandauMatrix(int k1,int k2,
					     int i1,int i2);
  virtual int whatTypeOfLandauMatrixDoINeed();
public:
  spinDownDensOperator(const Basis &newBasis,const eigSt &state_to_eval,int new_type);
  virtual ~spinDownDensOperator();
};

#endif
