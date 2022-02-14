/* Operator of the x-component of the total spin. In the present form
   it can only work with (end)basis yosBasis. */

#ifndef SPINXOPERATOR_HPP
#define SPINXOPERATOR_HPP


#include <Operator.hpp>
#include <yosBasis.hpp>

class spinXOperator : public Operator
{
private:
/* when computing <a|S_x|b>=\sum_{i,j} a_ib_j.<v_i|S_z|v_j>
  what at least a_ib_j has to be so that the matrix
  element of S_z would be computed (the matrix els' are 
  always of the order of unity. */
  const static double MIN_CONTRIB;   


/* endBases (in addition to Operator::mybasis). 
 */
  yosBasis *endYosBasis; 

/* Data (parameters) acquired from endBasis.
 */
  int Ne;
  int Nm;
  
/* This variable is to be set by ctor and should be 0 for standard
   usage of spinXOperator. In future, this operator might be also employed 
   by a new class spinZOperator (which would be similar to
   spinDensityOperator) which does the testing of equality of spatial
   parts itself (see spinSquareOperator::matEl_nonzero): 
   this is the only case when it is reasonable to set doNotTestEqSpat=1.
   */
  int doNotTestEqSpat;

/* Normalize the n-dimensional state in st_in and put the
   result into st_out. Returns 1 if the state is the zero 
   vector, otherwise returns 0. It is actually not necessary that this
   method is private... */  
  int normalize(double *st_out,double *st_in,int n);

/* The actual calculation of the S_x matrix element: between i_st1-th and
   i_st2-th state of endYosBasis. The result is put to *ME and the return 
   value is 1 if the MatEl was non-zero and 0 if it was zero.
*/
  int spinXMatEl(int i_st1,int i_st2,double *ME);

protected:
  const yosState *st1, *st2;

public:
  /* Standard ctor: doNotTestEqSpat=0. */
  spinXOperator(Basis *new_basis,int new_type);

  /* ctor for usage in spinXDensOperator: doNotTestEqSpat=0. */
  spinXOperator(Basis *new_basis);

private: 
  /* The main (common) part of all ctors. */
  void ctorBody();

public:

/* Computes the expectation value of S_x for the i_st-th vector 
   of my_basis. Can deal with non-normalized vectors too. */
  double totSpinX(int i_st);

/* Computes the expectation value of S_x for the vector 
   given by coefficients *inp_st with respect to endYosBasis.
   Can deal with non-normalized vectors too. 
   For the sake of security, dimension of *inp_st has to be
   given (which must match the dimension of endYosBasis.
*/  
  double totSpinX(double *coef,int check_dim);

  virtual int matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a); 
};


#endif













