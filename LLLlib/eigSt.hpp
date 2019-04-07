/*!
  \file eigSt.hpp
  \brief Holds the declaration of class eigSt
*/


#ifndef EIGST_HPP
#define EIGST_HPP

#include <myaccessories.hpp>
#include <Basis.hpp>
#include <iostream>


/*!
  \class eigSt
  \brief A class for keeping info about a state, presumably after diagonalization.
  
  It is intended as a heap of all thinkable data relevant to one state,
  i.e. objects of this class should be used to store small number of 
  states (100, 1000...) and not e.g. all states of a basis (10000, 100000...).
*/

class eigSt {

protected:
  double *coef;
  int dim;         // Dimension (i.e. size of coef*)
  bool m_bComplex; // wether eigenstate has complex coefficients
  

public:
  eigSt();               //! Create 'an empty' state (no WF can be stored)

  /*!
    \fn eigSt(Basis *new_bas, bool bComplex) 
    \brief Create a full state (with place for WF)
    \param new_bas basis these coefficients refer to
    \param bComplex wether or not the coefficients are to be stored as complex numbers
  */
  eigSt(Basis *new_bas, bool bComplex = false); // Create a full state (with place for WF)
  //! copy-constructor
  eigSt(const eigSt & eig);   

  virtual ~eigSt();

  /*!
    \fn void reinitialize(Basis *new_bas, bool bComplex)
    \param new_bas Basis to which coefficients refer
    \param bComplex wether to save coefficients as complex numbers
    \brief Initializes eigSt i.e. after array allocation
  */
  void reinitialize(Basis *new_bas, bool bComplex = true);
 //! returns coef[i], i=0..dim-1 for real, i=0..2*dim-1 for complex
  double outCoef(int i) const;            

  //! sets coef[i] to C, i=0..dim-1 for real, i=0..2*dim-1 for complex
  void inCoef(int i, double C);
  //! copies coef[0..dim-1] to out
  void outAllCoef(double *out) const;       
  //! copies coef[0..2*dim-1] to out
  void outAllCoef(std::complex<double> *out) const;  
  //! copies in to coef[0..dim-1]
  void inAllCoef(double *in);       
//! copies in to coef[0..2*dim-1]
  void inAllCoef(std::complex<double> *in);  
//! Delete the internal coefficents
  void deleteALlCoefficents();

  /*!
    \fn const double *getCoefPtr() const
    \brief Returns a const pointer on the coefficients
    \return Const ptr on coef
  */
  const double *getCoefPtr() const { return coef; };

  /*!
    \fn bool isComplex() const
    \brief Returns wether state contains complex coefficients.
  */
  bool isComplex() const { return m_bComplex; };
 
   //! gets energy
  double getEn () const;         
  //! sets energy
  void setEn (double En);        
  
  double getSz() const { return Sz; };


  double getTotS() const { return totS; };
  double getTotJ() const { return totJ; };
  int dimension() const;

  /*!
    \fn double normsq()
    \brief Calculates the norm^2 of this state.
    \return norm^2 of state
  */
  double normsq() const;

 

  /*!
    \fn double & operator[](int index)
    \param index Index of the element
    \brief Returns value of i-th entry in the vector.
    If eigSt is complex, two coefficients are composed to one complex number.
    i = 0..dim-1 in both cases, real or complex

    Throws an CxOutOfBoundsError for bad values of i
  */
  std::complex<double> operator[](int index);


  /*!
    \fn void setC_i(int i, std::complex<double> C)
    \brief Sets Coefficient C(i) of this vector.
    i = 0...dim-1. 
    \param C: Value of coefficient to set.
    \param i the index of the coefficent
    If Vector is real, only real part of C is mandatory.
  */    
  void setC_i(int i, std::complex<double> C);

  /*!
    \fn bool  operator= (const eigSt & eig)
    \brief assignment operator
  */
  eigSt & operator= (const eigSt & eig);

  /*!
    \fn bool operator== (const eigSt & eig)
    \brief compares two states to be _identical_
    needed by setupLandauMatrixForState
  */
  bool operator== (const eigSt & eig);

public:
    // Data
  double En;             //! Energy
  double Sz;             //! Total spin z-component
  double totJ;           //! Total momentum
  double totS;           //! Total spin
  


  
  //////////////////////
 // friend operators //
//////////////////////
 
  /*!
    \fn friend bool operator< (const eigSt &left, const eigSt &right)
    \brief Compares eigSt by energy (needed by list for sorting)
  */
  friend bool operator < (const eigSt &left, const eigSt &right);
};



#endif
