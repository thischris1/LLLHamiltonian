/*!
\file BarrYEl.cpp
*/

/*!
  \brief implementation of class BarrYEl
  stores matrix-elements of barrier parallel to x-axes (BarrY)
  and provides adressing of elements by means of operator(i,j)
*/

#include <BarrYEl.hpp>

// Constructor
// Nm = # of states 
BarrYEl::BarrYEl(int Nm) 
{
  m_Nm = Nm;
  m_pEl = new std::complex<double> [Nm];
  assert (m_pEl != NULL);
}

BarrYEl::~BarrYEl() 
{
  delete[] m_pEl;
}
 
/*
// general case of hermitian operator

BarrYEl::BarrYEl(int Nm) 
{
  m_Nm = Nm;
  m_pEl = new std::complex<double> [Nm*(Nm+1)/2];
  assert (m_pEl != NULL);
}


// adressing an Element
std::complex<double>  BarrYEl::operator()(int i, int j) 
{
  bool bConj = false;

  if (i < j) 
    {
      int tmp = i; 
      i = j; j = tmp;
      bConj = true;
    }
       
  // assert( j<=i )
  
  int at = (i*(i+1))/2 + j;

  // return original or conjugate  
  return (bConj) ? conj ( m_pEl[at] ) : m_pEl[at];
}

void BarrYEl::setAt(int i, int j, std::complex<double> val) 
{
  assert (i>=0 && i<m_Nm && j>=0 && j<m_Nm && j <= i);

  int at = (i*(i+1))/2 + j;
  m_pEl[at] = val;
}
*/

// special case only dependent of i-j
void BarrYEl::setAt(int i_j, std::complex<double> val) 
{
  assert ((i_j>=0) & (i_j<m_Nm));
  m_pEl[i_j] = val;
}

std::complex<double>  BarrYEl::operator()(int i, int j) 
{
  bool bConj = false;

  int i_j;

  if (i < j) 
    {
      i_j = j-i; 
      bConj = true;
    }
  else
  {
      
      i_j = i-j;
  }
  
  // assert( j<=i )
  
  // return original or conjugate  
  return (bConj) ? conj ( m_pEl[i_j] ) : m_pEl[i_j];
}
