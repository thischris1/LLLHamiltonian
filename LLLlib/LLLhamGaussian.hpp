/*! 
  \file LLLhamGaussian.hpp 
  \brief holds the declaration of class LLLhamGaussian
  \author jakob sachs und christian mueller
*/


#ifndef LLLhamGaussian_HPP
#define LLLhamGaussian_HPP

#include <iostream>
#include <Basis.hpp>
#include <yosBasis.hpp>
#include <Operator.hpp>
#include <LLLhamiltonian.hpp>
#include <BarrYEl.hpp>
#include <cassert>
#include <vector>
#include <geometry/CxImpurity.hpp>

typedef std::complex <double> dcomplex;

/*!
  \class LLLhamGaussian : public LLLhamiltonian 
  \brief storing matrix-elements forgaussian hole
  Class for representing a LLLhamiltonian with a gaussian shaped hole at a certain position.

*/
class LLLhamGaussian : public LLLhamiltonian
{
private:  
  // parameters controlling the barrier
 
    //!  energy-offset
  double m_dOffset;   

  //! matrix-el. for the hole (separated in two parts)
 
  double *m_pdHoleXeven;
  double *m_pdHoleXodd;

  BarrYEl *m_pHoleYeven;
  BarrYEl *m_pHoleYodd;
  int m_iBarrX;
   //! Impurities here
 //  std::vector<CxImpurity> impurities;
 
  int xKmax;
  int yKmax;
   CxImpurity m_impurity;
  ///////////////////////
 // auxiliary members //
///////////////////////

  


  /*!
  \fn   void calc_Hole (double s, double px, double py, double wx, double wy, 
  double a_to_b, double alpha1, double alpha2, int kxmax, int kymax);
  \brief Calculates matrix-elements for a gaussian-shaped hole
  \param alpha1 flux of solenoid 1 in units of h/e
  \param alpha2 flux of solenoid 2 in units of h/e
  */
  void calc_Hole ( double alpha1, double alpha2);

  /*!
    complex <double> calcMEXeven(double alpha1, double alpha2)
  */
  std::complex <double> calcMEXeven(double alpha1, double alpha2);

  /*!
    complex <double> calcMEXodd(double alpha1, double alpha2)
  */
  std::complex <double> calcMEXodd(double alpha1, double alpha2);
/*!
    complex <double> calcMEYeven(double alpha1, double alpha2)
  */
  std::complex <double> calcMEYeven(double alpha1, double alpha2);
  
  /*!
    complex <double> calcMEXeven(double alpha1, double alpha2)
  */
  std::complex <double> calcMEYodd(double alpha1, double alpha2);

  double iplusj(int i_plus_j,
		 double  alpha2,
		 double & xEvenPass, 
		 double & yEvenPass) const;


  bool iminusj(  std::complex <double> & MEYEvenPass,
		 std::complex <double> & MEYOddPass,
		 double  alpha1,
		 int  i_minus_j) const ;
   /*!
    \brief calculates the matrix elements of the Impurities on demand (no precalculation anymore)!
   */
  std::complex <double> impurityMatEl (int , int);


  /*!
    \brief calculcates the differnce dependent part of the matrix element
   */
  bool diffMatEl(int, std::complex <double>&, std::complex <double>&) const;

  /*!
    \brief calculcates the difference dependent part of the matrix element
  */
  bool sumMatEl(int, double &, double &) const;
  
  
  ////////////////////
 // public members //
////////////////////

public:
  /*!
    \fn LLLhamGaussian( yosBasis *basis, int mattype, double bli, double a, double b, double Strength, double xWidth, double xPosition, double holeStrength, double holeWidth, double holePosition, double offset, int iPotType, int xkmax, int yKmax )
    \brief Constructor fully qualified
    \param basis  used basis to calculate in
    \param type  type of matrices beeing used (full, sparse, complex) see Operator.hpp
    \param k_max  limit for the (infinite) sum
    \param strength  peak energy of potential in units of e^2/(epsilon*l0^2)  (l0 = magn. length)
    \param xWidth x- width of gaussian impurity in units of a
    \param xPosition  Position of the impurity in units of a.
    \param yWidth y-width of gaussian impurity in units of a
    \param yPosition Position of the impurity in units of a.
    \param offset  energy-offset (needed) for keeping 
    \param xKmax maximum k to sum over in matrix-elements' x-sums
    \param yKmax maximum k to sum over in matrix-elements' y-sums
  */
  LLLhamGaussian(yosBasis *basis, 
		 int type,
		 double bli,
		 double a,
		 double b,
		 double holeStrength,
		 double xWidth, 
		 double xPosition,
		 double yWidth,
		 double yPosition,
		 double offset,
		 int xkmax = 100,
		 int yKmax = 100);


LLLhamGaussian(yosBasis *basis,  int type, double bli,
	       double a, double b, const CxImpurity &  newImpurities,
	        double offset, int xkmax =100, int yKmax=100 );

  virtual ~LLLhamGaussian();
  
  /*!
    \fn virtual int matEl_perturb_nonzero(int i_st1, int i_st2, std::complex<double> *matEl)
    \brief Calculates matrix element between many-particle-states i_st1 and i_st2
    \param i_st1  index of many-particle state1 in basis passed in constructor
    \param i_st2  index of many-particle state2 in basis passed in constructor
    \param matEl  variable where value is returned
  */
  virtual int matEl_perturb_nonzero(int i_st1, int i_st2, std::complex<double> *matEl);

  /*!
    \fn addImpurity()

  */


  bool setImpurity(const CxImpurity &);


	CxImpurity getImpurity()const { return m_impurity;}
  /*!
    \brief calculates the contribution of impIndex impurity (matrix element iIndex, jIndex
    \param impIndex 
  */
  std::complex <double> singleME ( int  iIndex,
				   int  jIndex,
				   double alpha1,
				   double alpha2) const;

  std::complex <double> totalMatEl(int iIndex,
				   int jIndex,
				   double alpha1,
				   double alpha2) const;
				   
				   /*!
    \fn void dumpElements ()
    \brief Prints matrixelements of the hole on the screen
  */
  void dumpElements ();
};


#endif
