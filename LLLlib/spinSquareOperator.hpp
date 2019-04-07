/* Operator of the total spin (square). In the present form
   it can only work with (end)basis yosBasis. */

#ifndef SPINSQUAREOPERATOR_HPP
#define SPINSQUAREOPERATOR_HPP


#include <Operator.hpp>
#include <yosBasis.hpp>

/*!
  \class spinSquareOperator : public Operator
*/


class spinSquareOperator : public Operator
{
private:
/* when computing <a|S^2|b>=\sum_{i,j} a_ib_j.<v_i|S^2|v_j>
  what at least a_ib_j has to be so that the matrix
  element of S^2 would be computed (the matrix els' are 
  always of the order of unity. */
  const static double MIN_CONTRIB;   
  const static double MIN_CONTRIB_SINGLE; /* similar, 
			 if |a_i|<.., skip all b_j's */
		 


/* endBases (in addition to Operator::mybasis). 
 */
  yosBasis *endYosBasis; 

/* Data (parameters) acquired from endBasis.
 */
  int Ne;
  int Nm;
  
/* This variable is to be set by ctor and should be 0 for standard
   usage of spinSquareOperator. This operator is also employed by
   spinDensityOperator which does the testing of equality of spatial
   parts itself (see spinSquareOperator::matEl_nonzero): 
   this is the only case when it is reasonable to set doNotTestEqSpat=1.
   */
  int doNotTestEqSpat;

/* Normalize the n-dimensional state in st_in and put the
   result into st_out. Returns 1 if the state is the zero 
   vector, otherwise returns 0. It is actually not necessary that this
   method is private... */  
  int normalize(double *st_out,double *st_in,int n);

/* The actual calculation of the S^2 matrix element: between i_st1-th and
   i_st2-th state of endYosBasis. The result is put to *ME and the return 
   value is 1 if the MatEl was non-zero and 0 if it was zero.
*/
  int spinMatEl(int i_st1,int i_st2,double *ME);

protected:
  const yosState *st1, *st2;

public:
  /* Standard ctor: doNotTestEqSpat=0. */
  spinSquareOperator(Basis *new_basis,int new_type);

  /* Ctor for the case when doNotTestEqSpat=1. */
  spinSquareOperator(Basis *new_basis);
private: 
  /* The main (common) part of all ctors. */
  void ctorBody();

public:

/* Computes the expectation value of S^2 for the i_st-th vector 
   of my_basis. Can deal with non-normalized vectors too. */
  double totSpin(int i_st);

/* Computes the expectation value of S^2 for the vector 
   given by coefficients *inp_st with respect to endYosBasis.
   Can deal with non-normalized vectors too. 
   For the sake of security, dimension of *inp_st has to be
   given (which must match the dimension of endYosBasis.
*/  
  double totSpin(double *coef,int check_dim);

  virtual int matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a); 
};


#endif













