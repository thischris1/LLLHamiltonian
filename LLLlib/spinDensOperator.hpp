/*!
  \file spinDensOperator.hpp
*/
#ifndef SPINDENSOPERATOR_HPP
#define SPINDENSOPERATOR_HPP

#include <densOperator.hpp>
#include <spinSquareOperator.hpp>
#include <fstream>
#include <iostream>

class spinDensOperator : public DensOperator
{
private:
  spinSquareOperator *spSqOp;
protected:
  virtual std::complex<double> getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2);
  virtual int whatTypeOfLandauMatrixDoINeed();
public:
  spinDensOperator(const Basis &newBasis,const eigSt &state_to_eval,int new_type);
  virtual ~spinDensOperator();
  virtual int s1_nr(int i) const;
  virtual int s2_nr(int i) const;
};

#endif
