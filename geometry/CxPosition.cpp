/*!
\file CxPosition.cpp
  \brief Contains implementation of class CxPosition
  \author Christian Mueller
  \date 07 May 03

*/
#include <geometry/CxPosition.hpp>
#include <utils/logger.hpp>

/*!
  \brief standard ctor postion is 0.0 afterwards
	
*/

CxPosition::CxPosition():
xPos(0.0),
yPos(0.0),
epsilon(1e-3)
{
  glLogger.debug("Standard ctor of CxPosition::CxPosition");
}


/*!
  \param newX x Position
  \param newY y-Position
  \brief ctor from coordinate values
*/
CxPosition::CxPosition(double newX, double newY)
:xPos(newX),
yPos(newY),
epsilon(1e-10)
{
    setXPosition(newX);
    setYPosition(newY);
}

/*!
  \return the x position 

*/
double CxPosition::getXPosition(void) const
{ 
  return (xPos);
}
/*!
  \return the y-position

*/
double CxPosition::getYPosition(void) const
{
  return (yPos);
}

/*!
  \brief Set X Position ( ignored if < 0.0)
  \param newXpos x-Position 
*/
void CxPosition::setXPosition (double newXpos)
{ 

      xPos = newXpos;

}
/*!
  \brief Set Y Position ( ignored if < 0.0)
  \param newYpos the new y-position value
  \todo Think about vectors < 0
*/

void CxPosition::setYPosition (double newYpos)
{
   yPos = newYpos;
}


/*!
  \return true if positions do match (within epsilon)
*/

bool CxPosition::compare(const CxPosition &rhs) const
{
  glLogger.debug("Comparing x1 =(%f), x2 = (%f), y1 = (%f), y2 = (%f) with epsilon (%f)",
		 xPos,  rhs.getXPosition(), yPos, rhs.getYPosition(), epsilon);
    if (abs(xPos - rhs.getXPosition()) <= epsilon)
    {
      if (abs(yPos - rhs.getYPosition()) <= epsilon)
	{
	  glLogger.debug("Same positions");
	  return (true);
	}
    }
    glLogger.debug("Different positions");
  return false;

}

/*!
  \param newValue the new epsilon (should be > 0)
*/
void CxPosition::setEpsilon(double newValue )
{
  if (newValue > 0)
    {
      epsilon = newValue;
    }

}

/*!
  \return the current Epsilon value
*/

double CxPosition::getEpsilon(void) const
{

  return (epsilon);


}

/*!
\brief convenience operator
\return true if equal (within epsilon)
\param rhs the other Position to be comapred with
*/

bool CxPosition::operator == (const CxPosition &rhs)
{


  return (compare(rhs));
}

/*!
  \brief Difference vector a - b
  \param rhs the vector b
  \return the difference vector
*/
const CxPosition  CxPosition::operator -(const CxPosition &rhs)
{
  glLogger.debug("Difference vector between (%f), (%f) and (%f) (%f)",
		xPos, yPos,  rhs.getXPosition(), rhs.getYPosition());
  double newX = xPos - rhs.getXPosition();
  double newY = yPos - rhs.getYPosition();
  glLogger.info("Differenz (%f), y(%f)", newX, newY);
  CxPosition retVal(newX, newY);
  glLogger.debug("Result is (%f), (%f)", retVal.getXPosition(), retVal.getYPosition());
  return(CxPosition(retVal.getXPosition(), retVal.getYPosition()));

}

/*!
  \return the euclidean norm of the vector
*/
double CxPosition::vabs(void)const
{

  double retVal= 0;
  retVal = xPos*xPos + (yPos*yPos);
  retVal = sqrt(retVal);
  return (retVal);

}


/*!
	\fn double  CxPosition::vabs(CxPosition inValue)
	\author Christian Mueller
	\date  27 Aug 03
	\return 
	\param CxPosition inValue
 	\brief Convenience function  
 */
/*
	Pre	: 
	Post	: 

 */ 

double CxPosition::vabs(CxPosition inValue) const
{
    // Test of Parameter

    return (inValue.vabs());
    

}

/*!
\fn void CxPosition::multiplyVector(double newFac)
\brief vector times scalar multiplication. Result changes the current vector
\param newFac factor by which the vector shall be enlarged, shortened

*/
void CxPosition::multiplyVector(double newFac)
{
  xPos = xPos*newFac;
  yPos = newFac*yPos;

}


/*!
	\fn CxPosition  CxPosition::operator+(CxPosition &rhs)
	\author Christian Mueller
	\date  26 Aug 03
	\return 
	\param CxPosition &rhs
 	\brief  Vector addition
*/

const CxPosition  CxPosition::operator +(const CxPosition &rhs) 
{
    
  return ( CxPosition ((xPos + rhs.getXPosition()),
		       (yPos + rhs.getYPosition())));


}

/*!

	\author Christian Mueller
	\date  27 Aug 03
	\return index of the nearest Position to this Position  
	\param vPositions  vector of CxPositions to be searched through
 	\brief Searches the nearest Position to the current object.

        May be used with electronSites to find the electron a vortex belongs to.
        
 */


int CxPosition::getNearestPosition(std::vector<CxPosition> vPositions)
{
    // Test of Parameter
  glLogger.info("Entering getNearestPosition wiuth positions (%f), (%f)",
		xPos, yPos);
    int retVal = 0;
    double minimalDistance = 100; // There must be a start value
    
    if (vPositions.size() == 0)
    {

        return (retVal);
    }
    
    for (unsigned int vIndex = 0; vIndex < vPositions.size(); vIndex++)
    {
      CxPosition differenzVec = (*this) - vPositions[vIndex];
      
      double tempDistance =  vabs((*this)-vPositions[vIndex]);
      glLogger.info("Looking up Position with index (%d), at position (%f), (%f)",
		    vIndex,
		    vPositions[vIndex].getXPosition(),
		    vPositions[vIndex].getYPosition());
      glLogger.info("DifferenzVector is  (%f), (%f), Betrag is (%f)",
		    differenzVec.getXPosition(),
		    differenzVec.getYPosition(),
		    differenzVec.vabs());
      glLogger.debug("Abstand is (%f)", tempDistance);
      if (tempDistance < minimalDistance)
        {
	  retVal = vIndex;
	  minimalDistance =  tempDistance;
        }
      
    }
    return (retVal);
    

}


/** No descriptions */
 CxPosition::~CxPosition(){
}
/** Copy ctor */
 CxPosition::CxPosition(const CxPosition &rhs)
: xPos(rhs.getXPosition()),
yPos(rhs.getYPosition()),
epsilon(rhs.getEpsilon())
{
 
}

 CxPosition & CxPosition::operator =(const CxPosition &rhs)
{
  if (this != &rhs)
    {
      xPos = rhs.getXPosition();
      yPos = rhs.getYPosition();
    }
  return (*this);
}
