/*!
\file BarrYEl.hpp
*/

/*!
  \class BarrYEl
  \brief implementation of class BarrYEl
  stores matrix-elements of barrier parallel to x-axes (BarrY)
  and provides adressing of elements by means of operator(i,j)
*/

#ifndef BARRYEL_HPP
#define BARRYEL_HPP

#include <complex>
#include <cassert>

class BarrYEl 
{

private:
  int m_Nm;
  std::complex<double> *m_pEl;

public:

  /*!
  \brief constructor
  Nm = number of states
  */
  BarrYEl(int Nm);

  /*!
    \brief destructor
  */
  ~BarrYEl();

  /*!
    adressing of Element <i|V|j>
  */
  std::complex<double> operator()(int i, int j);
  
  /*!
    \brief setting elements
    i_j = i-j, because matrix elements only depend on i-j here.
    assumed that i>=j, i-j in [0,Nm-1]
  */
  void setAt (int i_j, std::complex<double> val);

//void setAt (int i, int j, std::complex<double> val);
};

#endif
