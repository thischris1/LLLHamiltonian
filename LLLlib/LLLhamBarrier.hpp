/*! 
  \file LLLhamBarrier.hpp 
*/


#ifndef LLLhamBarrier_HPP
#define LLLhamBarrier_HPP

#include <iostream>
#include <Basis.hpp>
#include <yosBasis.hpp>
#include <Operator.hpp>
#include <LLLhamiltonian.hpp>
#include <BarrYEl.hpp>
#include <cassert>

/*!
  \class LLLhambBarrier : public LLLhamiltonian 
  \brief storing matrix-elements for barrier || x-axes
  Class for representing a LLLhamiltonian with a gaussian shaped barrier with a gaussian hole
  and delta shaped barrier with a gaussian hole
*/
class LLLhamBarrier : public LLLhamiltonian
{
private:  
  // parameters controlling the barrier
  int m_iPotType;       // Type of potential in x-direction: 0->gaussian, 1->delta
  double m_iBarrX;
  double m_dOffset;    // energy-offset
  double m_dBarrStrength;
  
  // matrix-el. calculated by calc_BarrierX_el for barrier in x-direction
  double *m_pdBarrXEl;

  // matrix-el. claculated by calc_BarrierY_el for barrier in y-direction
  BarrYEl *m_pBarrYEl;

  // matrix-el. for the hole (separated in two parts)
  // only depends on 
  double *m_pdHoleXeven;
  double *m_pdHoleXodd;
  BarrYEl *m_pHoleYeven;
  BarrYEl *m_pHoleYodd;

 
  ///////////////////////
 // auxiliary members //
///////////////////////

  /*!
    \fn   void calc_BarrierX_el (double p, double s, double w, double alpha2, int kmax)
     \brief evalutates one-particle matrix-elements of barrier along x-axis, parralel to y-axis
     \param p  Position of barrier in units of a
     \param s  peak strength in units of e^2/(epsilon l0)
     \param w  width in units of a
     \param alpha2 flux of solenoid 2 in units of h/e modifying periodic boundary conditions
     \param kmax maximum k to sum over
  */
  void calc_BarrierX_el (double p, double s, double w, double alpha2, int kmax);

  /*!
    \fn   void calc_BarrierY_el (double p, double s, double w, double alpha1, int kmax)
     \brief Evalutates one-particle matrix-elements of barrier along y-axis, parrallel to x-axis
     \param p  Position in units of b
     \param s  peak strength in units of e^2/(epsilon l0)
     \param w  width in units of b
     \param alpha1 flux of solenoid1 affecting x-momentum
     \param kmax maximum k to sum over
  */
  void calc_BarrierY_el (double p, double s, double w, double alpha1, int kmax);

  /*!
  \fn   void calc_Hole (double s, double px, double py, double wx, double wy, 
  double a_to_b, double alpha1, double alpha2, int kxmax, int kymax);
  \brief Calculates matrix-elements for a gaussian-shaped hole
  \param s  peak strength in units of e^2/(epsilon l0)
  \param px  x-position in units of a
  \param py  y-position in units of b
  \param wx  width of gaussian in x-direction in units of a (=2*sigma)
  \param wy  width of gaussian in y-direction in units of b (=2*sigma)
  \param a_to_b  aspect ratio
  \param alpha1 flux of solenoid 1 in units of h/e
  \param alpha2 flux of solenoid 2 in units of h/e
  \param kxmax   maximum k to sum over in x-sum
  \param kymax   maximum k to sum over in y-sum
  */
  void calc_Hole ( double s, 
		   double px, double py, 
		   double wx, double wy, 
		   double a_to_b, 
		   double alpha1, double alpha2,
		   int kxmax, int kymax);

  /*!
    \fn void dumpElements ()
    \brief Prints matrixelements of the hole on the screen
  */
  void dumpElements ();

  ////////////////////
 // public members //
////////////////////

public:
  /*!
    \fn LLLhamBarrier( yosBasis *basis, int mattype, double bli, double a, double b, double Strength, 
    double xWidth, double xPosition, double holeStrength, double holeWidth, double holePosition, double offset, 
    int iPotType, int xkmax, int yKmax )
    \brief Constructor
    \param basis  used basis to calculate in
    \param mattype  type of matrices beeing used (full, sparse, complex) see Operator.hpp
    \param k_max  limit for the (infinite) sum
    \param strength  peak energy of potential in units of e^2/(epsilon*l0^2)  (l0 = magn. length)
    \param width  width of gaussian profile in units of a
    \param xposition  Position of the barrier in units of a. If iPotType = 1 the barrier is located on the next Xj
    \param offset  energy-offset (needed) for keeping 
    \param alpha1 flux of solenoid1 affecting x-momentum in units of h/e
    \param alpha2 flux of solenoid2 affecting y-momentum in units of h/e
    \param iPotType  0 -> gaussian shaped barrier; 1 -> delta shaped barrier
    \param xKmax maximum k to sum over in matrix-elements' x-sums
    \param yKmax maximum k to sum over in matrix-elements' y-sums
  */
  LLLhamBarrier(yosBasis *basis, 
		int type,
		double bli,
		double a,
		double b,
		double Strength, 
		double xWidth, 
		double xPosition,
		double holeStrength,
		double holeWidth,
		double holePosition,
		double offset,
		double alpha1,
		double alpha2,
		int iPotType,
		int xkmax,
		int yKmax );

  ~LLLhamBarrier();
  
  /*!
    \fn virtual int matEl_perturb_nonzero(int i_st1, int i_st2, std::complex<double> *matEl)
    \brief Calculates matrix element between many-particle-states i_st1 and i_st2
    \param i_st1  index of many-particle state1 in basis passed in constructor
    \param i_st2  index of many-particle state2 in basis passed in constructor
    \param matEl  variable where value is returned
  */
  virtual int matEl_perturb_nonzero(int i_st1, int i_st2, std::complex<double> *matEl);

};


#endif
