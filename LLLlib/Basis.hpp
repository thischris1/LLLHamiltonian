/*!
  \file Basis.hpp
  \brief Declaration of class Basis
  \author Karel Vyborny
*/
#ifndef BASIS_HPP
#define BASIS_HPP

#include <State.hpp>
#include <stdio.h>
#include <myaccessories.hpp>
#include <ERRORS.h>
#include<complex>
#include <vector>
typedef enum {LINEAR_COMBINATION,YOS_BASIS} BASIS_TYPE;
/*!
  \class Basis
  \brief abstarct Baseclass for many particle bases. Do not instantiate it directly but use derived classes.
  \author Karel Vyborny
*/
class Basis 
{


public:
  /*!
    \fn Basis(Basis *new_basis,double *coef,int n)
    \brief    Initialize a basis (or better) a set of n vectors 
     defined by coef with respect to some other basis 
     (new_basis).
     \param n 
     \param coef The coef array is sorted as n clusters (vectors) of
     the length new_basis->dimension().
  */
  Basis(Basis *new_basis,double *coef,int n);

// copy ctor
	Basis(const Basis &rhs);
  //! Access to the i-th state of the basis. 
  virtual State *operator[](long int i) const =0;

  //! Guess what.
  long int dimension() const;

//!
	void setDimension(int n_dim);
  //! Get the i-th component of the n_state-th vector.
  virtual double getCoef(int n_state,int i);

  /*!

  \fn  virtual double norm(int i_st)
  \brief Computes the norm of the i-th vector;
  CAUTION: do not forget to redefine this method in 
  derived classes (e.g. a yosState, i.e. a product
  state of one-electron normed WF's has automatically
     norm=1. 
  */
  virtual double norm(int i_st);

  //! With respect to which basis is the new one 
  virtual Basis *getMyBasis()const;
  //!  Destructor
  virtual ~Basis();
  
  //Basis operator = (const Basis &rhs);
  BASIS_TYPE getBasisType() const;
  
  /* An object of Basis type is a basis described as a linear
     combination of vectors of some other basis which are LC
     of yet another basis etc. This sequence should end
     by an 'implemented' basis (e.g. yosBasis). Get a pointer
     to this basis. */
  virtual Basis *getEndBasis();

  /* Value of SP wavefct at place (x,y). Not thought to be
     implemented here (gives an err_msg). It is implemented
     e.g. in yosBasis. 
  virtual complex<double> SPwaveFct(int i,double x,double y);
  */
  	double *getAllCoefs () const;

	std::vector <double> getAllVCoefs() const;

protected:
  long int dim;
  Basis *myBasis;
  BASIS_TYPE basis_type;
 
 private:	
  double *coef;
  std::vector<double> m_coef;
};


#endif









