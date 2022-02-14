/*!
  \file <persist_eigSt.cpp
  
*/
#include <persist_eigSt.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>


//persist_eigSt::persist_eigSt() : eigSt () {
//}
/*!
  \fn persist_eigSt::persist_eigSt(Basis *pBasis_to_use, bool bComplex)
  \param pBasis_to_use pointer to Basis
  \param bComplex true if coefficents are complex
  \brief Constructtor

*/

persist_eigSt::persist_eigSt(Basis *pBasis_to_use, bool bComplex)
: eigSt (pBasis_to_use, bComplex),
  m_bBinary(false)

{
}


/*!
  \fn persist_eigSt::persist_eigSt (const persist_eigSt & src)
  \param src to be copied
  \brief Copy constructor

*/
persist_eigSt::persist_eigSt (const persist_eigSt & src) : eigSt (src)
{  
  m_bBinary = src.is_binary();
  

}

/*!
  \fn persist_eigSt::~persist_eigSt() 
  \brief Destructor, does nothing
*/
persist_eigSt::~persist_eigSt() 
{
  glLogger.info("persist eigSt destructor");
}

void persist_eigSt::set_binary(bool bYes) 
{
  m_bBinary = bYes;
}


/*
void persist_eigSt::operator= (const eigSt &eig) {
  
}
*/
/*!
  \fn std::ostream & operator<< (std::ostream &out, const persist_eigSt &eigst) \brief Related fucntion to allow streaming of coefficents 
*/
std::ostream & operator<< (std::ostream &out, const persist_eigSt &eigst) 
{

  // if (eigst.isComplex())
  //glLogger.info ("Saving complex vector");

  int nCount = eigst.isComplex() ? 2*eigst.dimension() : eigst.dimension();

  if ( eigst.is_binary() )
    {
      // convert doubles to unformatted binary string
      double d = eigst.getEn();
      out.write( (char*) &d, sizeof(double)/sizeof(char) );
      out.write ( (char*) eigst.getCoefPtr(), nCount*sizeof(double)/sizeof(char) );      
    }
  else
    {
      out << eigst.getEn();
      for (int i=0; i< nCount; i++) 
	out << ' ' << eigst.outCoef(i);
      out << '\n';
    }

  return out;
}

/*!
  \fn std::istream & operator>> (std::istream &in, persist_eigSt &eigst) 
  \brief Streams the content of state to stream
*/
std::istream & operator>> (std::istream &in, persist_eigSt &eigst) 
{
  double En = 0.0;

  if (eigst.is_binary()) 
    {
      int nCount = eigst.isComplex() ? 2*eigst.dimension() : eigst.dimension();
      double *pCoefs = new double [nCount];
      


      in.read( (char*) &En, sizeof(double)/sizeof(char) );
      eigst.setEn(En);
      
      in.read ( (char*) pCoefs, sizeof(double)*nCount/sizeof(char) );
      
      if (eigst.isComplex())
	eigst.inAllCoef( (std::complex<double>*)pCoefs );
      else
	eigst.inAllCoef(pCoefs);
	 
      delete [] pCoefs;
    }
  else
    {
      in >> En;
      eigst.setEn(En);

      int nCount = eigst.isComplex() ? 2*eigst.dimension() : eigst.dimension();
      double dCoefi;
      for (int i = 0; i < nCount; i++) 
	{
	  in >> dCoefi;
	  eigst.inCoef(i, dCoefi);
	}
    }
  
  return in;
}



