/*!
  \file eigSt.cpp
  \brief Holds implementation of class eigSt, which is a class for keeping info about a state, presumably after diagonalization 
  \author Karel Vyborny, Christian Mueller

*/

#include <eigSt.hpp>
#include <cassert>
#include <utils/logger.hpp>
#include <ERRORS.h>




eigSt::eigSt(Basis *new_bas, bool bComplex):
coef(0),
dim(new_bas->dimension()),
m_bComplex(bComplex),
En(0),
Sz(0),
totJ(0),
totS(0)
{
 
   try {
    if (bComplex)
      {
	coef= new double [2*dim];
      }
    else
      {
	coef = new double [dim];
      }
  }
  catch (...)
    {
      coef = 0;
      throw (CxErrors("Out of memory in constructor of eigSt"));
    }



}


/*!
  \fn eigSt::eigSt(const eigSt & eig)

*/
eigSt::eigSt(const eigSt & eig)
{
  glLogger.info("copy constructor eigSt");

  m_bComplex = eig.isComplex();

  En = eig.getEn();
  Sz = eig.getSz();
  totS = eig.getTotS();
  totJ = eig.getTotJ();

  dim = eig.dimension();

  if (m_bComplex)
    {
      coef = new double [2*dim];
      eig.outAllCoef ((std::complex<double>*)coef);
    }
  else
    {
      coef = new double[dim];
      eig.outAllCoef (coef);
    }

  //for (int i = 0; i < dim; i++)
  //  coef[i] = eig.outCoef(i);
  //
}

eigSt & eigSt::operator= (const eigSt & eig)
{
  //std::cout << "assignment operator eigSt\n";
  if (this != &eig)
    {
      
      En = eig.getEn();
      Sz = eig.getSz();
      totS = eig.getTotS();
      totJ = eig.getTotJ();
      dim = eig.dimension();
      m_bComplex = eig.isComplex();
      
      if (coef != NULL) delete [] coef;
      
      if (m_bComplex)
	{
	  coef = new double [2*dim];
	  eig.outAllCoef ((std::complex<double>*)coef);
	}
      else
	{
	  coef = new  double [dim];
	  eig.outAllCoef (coef);
	}

  //for (int i = 0; i < dim; i++)
  //  coef[i] = eig.outCoef(i);
    }

  return *this;
}

eigSt::eigSt():
coef(0),
dim(0),
m_bComplex(false),
En(0),
Sz(0),
totJ(0),
totS(0)
{
  
}

void eigSt::reinitialize(Basis *new_bas, bool bComplex )
{ 
  if (coef != NULL) delete[] coef;
  dim = new_bas->dimension();

  m_bComplex = bComplex;

  if (m_bComplex)
    coef = new double [2*dim];
  else
    coef = new double[dim];
}

eigSt::~eigSt()
{
  if (coef != NULL) 
    delete[] coef;
}

double eigSt::outCoef(int i) const
{
  int up_bound = m_bComplex ? 2*dim : dim;
  if(i<0 || i>=up_bound) 
    {
      throw CxOutOfBoundsError(__FILE__, __LINE__);
    };
  return coef[i];
}

void eigSt::inCoef(int i, double C)
{
  int up_bound = m_bComplex ? 2*dim : dim;
  if(i<0 || i>=up_bound) 
    {
      throw CxOutOfBoundsError(__FILE__, __LINE__);
    }
  coef[i] = C;
}



void eigSt::outAllCoef(double *out) const
{
  // if this assertion fails, you tried to read out the coefficients
  // of a complex vector in a double array
  if (m_bComplex)
    {
      out = 0;
      throw CxErrors ("you tried to read out the coefficients of a complex vector in a double array");
    }
  //int nCount = m_bComplex ? 2*dim : dim;
  for(int i=0; i<dim; i++)
    out[i] = coef[i];
}

void eigSt::outAllCoef(std::complex<double> *out) const
{
  // if this assertion fails, you tried to read out the coefficients
  // of a real vector in a complex array
  if (!m_bComplex)
  {
      out = 0;
      throw (CxErrors("you tried to read out the coefficients  of a real vector in a complex array"));
  }
  for(int i=0; i<dim; i++)
    out[i] = std::complex<double>(coef[2*i], coef[2*i+1]);
}

/*!
  \fn void eigSt::inAllCoef(double *in)
*/
void eigSt::inAllCoef(double *in)
{
 
   if (in == 0)
    {
      coef = 0;
      throw CxNullPointerError(__FILE__, __LINE__);
    }
    glLogger.info("eigSt read in inAllCoef(double)");
  for(int i=0; i<dim; i++)
    coef[i] = in[i];
}

void eigSt::inAllCoef(std::complex<double> *in)
{
  // if this assertion fails, you tried to read in the coefficients
  // of a real vector from a complex array
  if (in == 0)
    {
      coef = 0;
      throw CxNullPointerError(__FILE__, __LINE__);
    }
  
   glLogger.info("eigSt read in inAllCoef(double)");
  for(int i=0; i<dim; i++)
    {
    	if(coef)
    	{
    		/*
    		 * Delete old and reallocate
    		 */
    	}
    	
      coef[2*i] = in[i].real();
      coef[2*i+1] = in[i].imag();
    }
}

void eigSt::setEn (double new_En) {
  En = new_En;
}

double eigSt::getEn() const
{
  return En;
}


int eigSt::dimension() const 
{
  return dim;
}

double eigSt::normsq() const
{
  double normsq = 0.0;
  if (m_bComplex)    
    for (int i=0; i<dim; i++)
      {
	std::complex<double> z = ((std::complex<double>*) coef)[i];
	normsq += (conj(z)*z).real();
      }
  else
    for (int i=0; i<dim; i++)
      {	
	normsq += coef[i]*coef[i];
      }

  return normsq;
}


  ///////////////
 // operators //
///////////////

bool operator< (const eigSt &left, const eigSt &right)
{
  return (left.En < right.En);
}

/*!
  \fn std::complex<double> eigSt::operator[](int index)

*/
std::complex<double> eigSt::operator[](int index)
{

  if (index < 0 )
    {
      throw CxIndextoSmallError(index);
    }
  if (index > dim)
    {
      throw CxIndextoLargeError(index);
    }

  if (m_bComplex)
    return ( (std::complex<double>*)coef ) [index];
  else
    return std::complex<double> (coef[index], 0.0);

}
/*!
  \fn void eigSt::setC_i(int i, std::complex<double> C)
*/

void eigSt::setC_i(int index, std::complex<double> C)
{
   if (index < 0 )
    {
      throw CxIndextoSmallError(index);
    }
  if (index > dim)
    {
      throw CxIndextoLargeError(index);
    }

  if (m_bComplex)
    ( (std::complex<double>*)coef ) [index] = C;
  else
    coef[index] = C.real();
}


bool eigSt::operator== (const eigSt & eig)
{

  bool bEqual = 
    (m_bComplex == eig.isComplex())
    && (dim == eig.dimension())
    && (En == eig.getEn()) 
    && (totJ == eig.getTotJ()) 
    && (totS == eig.getTotS())
    && (Sz == eig.getSz());
  
  if (bEqual)
    {
      int nCount = m_bComplex ? 2*dim : dim;
      for (int i=0; i<nCount; i++)
	bEqual &= (coef[i] == eig.outCoef(i));
    }

  return bEqual;
}

void eigSt::deleteALlCoefficents()
{
	if (!coef)
	{
		/*
		 * Nothing to do
		 */
		return;	
	}	
	delete [] coef;
	if (m_bComplex)
	{
    	coef = new double [2*dim];
	}
  	else
  	{
    	coef = new double[dim];
  	}
	return;
}

