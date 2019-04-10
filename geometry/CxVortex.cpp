#include <geometry/CxVortex.hpp>
#include <LLLlib/LLLlib.h>
/*!


*/
CxVortex::CxVortex():
CxPeriodicPosition(),
windNumber(0)
{

}

/*!


*/
CxVortex::CxVortex(double xPos, double yPos, int windingNumber, double xSize, double ySize): CxPeriodicPosition(xPos,yPos, xSize, ySize)
{


  setWindingNumber(windingNumber);
}

CxVortex::CxVortex(CxPeriodicPosition pos, int windingNumber):CxPeriodicPosition(pos)
{
  setWindingNumber(windingNumber);

}


/*!
\brief Accessor
\return Winding number as multiple of 2\pi
*/
int CxVortex::getWindingNumber(void) const
{
  return (windNumber);

}

/*!
\brief setter of windingnumber
*/

void CxVortex::setWindingNumber (int newVal)
{
  /*
if (newVal < 0)
    {
      throw (CxBadValueError(__FILE__,__LINE__, "windingNumber is less then zero") );
      return;
    }
  */
 windNumber = newVal; 
}
/** No descriptions */
CxVortex CxVortex::operator =(const CxVortex & rhs){
	if (this != &rhs)
	{
		CxPosition::operator =(rhs);
		windNumber = rhs.getWindingNumber();

	}
	return (*this);
}
CxVortex::CxVortex(const CxVortex & rhs):
CxPeriodicPosition(rhs)
{
  windNumber = rhs.getWindingNumber();

}

const CxVortex CxVortex::operator-( const CxVortex & rhs)
{
  if (rhs.getWindingNumber() != windNumber)
    {
      throw (CxBadValueError(__FILE__,__LINE__, "windingNumber is not the same") );
    }
  return ( CxVortex(CxPeriodicPosition::operator-(rhs), windNumber));
  
}
