/***************************************************************************
                          cxperiodicposition.cpp  -  description
                             -------------------
    begin                : Sat Oct 23 2004
    copyright            : (C) 2004 by Christian Mueller
    email                : cmueller@physnet.uni-hamburg.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <geometry/CxPeriodicPosition.hpp>
#include <LLLlib/LLLlib.h>
/** Default ctor, private does nothing*/
CxPeriodicPosition::CxPeriodicPosition():
CxPosition(),
m_xCellSize(1.0),
m_yCellSize(1.0){
}

CxPeriodicPosition::CxPeriodicPosition(double  n_X, double  n_Y, double  aTob):
  CxPosition(xPos, yPos),
  m_xCellSize(1.0),
  m_yCellSize(aTob)
{
  xPos = mapXPos(n_X);
  yPos = mapYpos(n_Y);
  set_xCellSize(1.0);
  set_yCellSize(aTob);
}


/** Destructor */
CxPeriodicPosition::~CxPeriodicPosition(){
}
/** Read property of double m_xCellSize. */
double CxPeriodicPosition::get_xCellSize()const {

	return m_xCellSize;
}
/** Write property of double m_xCellSize. */
void CxPeriodicPosition::set_xCellSize( const double& _newVal){
	glLogger.debug("CxPeriodicPosition():setXCellSize set xSize to (%f)", _newVal);
	if (_newVal < 1e-06) {
	  throw CxBadValueError("x Cellsize 0 or negative", __LINE__, __FILE__);
	}
	m_xCellSize = _newVal;
}
/** Read property of double  m_yCellSize. */


CxPeriodicPosition::CxPeriodicPosition(const CxPeriodicPosition &rhs):
CxPosition(rhs)
{
  set_xCellSize(rhs.get_xCellSize());
  set_yCellSize(rhs.get_yCellSize());

}
/** No descriptions */
 CxPeriodicPosition::CxPeriodicPosition(double n_X, double n_Y, double n_XSize, double n_ySize):
	CxPosition(n_X, n_Y)
 {
 set_xCellSize( n_XSize);
 set_yCellSize( n_ySize);
 xPos = mapXPos(n_X);
 yPos = mapYpos(n_Y);
 glLogger.debug("CxPeriodicPosition with x=(%f), y=(%f), cell = (%f)x(%f)", xPos, yPos, m_xCellSize, m_yCellSize);
}
/** Maps a x Position into the unit cell */
double CxPeriodicPosition::mapXPos(double xPos) const {

 double retVal;
  if ( (xPos > 0  ) && (xPos < m_xCellSize))
    {
      return (xPos);
    }
  retVal = xPos -((floor(xPos/m_xCellSize))*m_xCellSize);
	glLogger.debug("mapXPos maps (%f) to (%f)",xPos, retVal);
return (retVal);

}
/** Maps the y Position into unit cell */
double CxPeriodicPosition::mapYpos(double yPos) const{
 double retVal;
  if ( (yPos > 0  ) && (yPos < m_yCellSize))
    {
      return (yPos);
    }
  retVal = (yPos - (floor(yPos/m_yCellSize))*m_yCellSize);
  	glLogger.debug("mapYPos maps (%f) to (%f)",yPos, retVal);	
return (retVal);



}
/** Equal if cell and Position identical */
bool CxPeriodicPosition::operator == (const CxPeriodicPosition &rhs)
{
	if 	(rhs.get_xCellSize() == m_xCellSize && rhs.get_yCellSize() == m_yCellSize)
	{
	  return (CxPosition::operator==(rhs));

	}
	else
	{
		return (false);
	}
}
/** Read property of double  m_YCellSize. */
double CxPeriodicPosition::get_yCellSize() const{
	return m_yCellSize;
}
/** Write property of double  m_YCellSize. */
void CxPeriodicPosition::set_yCellSize( const double & _newVal){
glLogger.debug("CxPeriodicPosition():setYCellSize set ySize to (%f)", _newVal);
 	if (_newVal < 1e-06) {
	  throw CxBadValueError("y Cellsize 0 or negative", __LINE__, __FILE__);
	}
 
	m_yCellSize = _newVal;
}

/** No descriptions */
CxPosition CxPeriodicPosition::mapPosition(CxPosition & inPos){

	CxPosition retVal(mapXPos(inPos.getXPosition()), mapYpos(inPos.getYPosition()));
	return (retVal);
}
/** Addition operator */
const CxPeriodicPosition  CxPeriodicPosition::operator + (const CxPeriodicPosition  &rhs){
if ( (rhs.get_xCellSize() != m_xCellSize) || (rhs.get_yCellSize() != m_yCellSize))
{
	throw (CxErrors("Bad cell sizes in CxPeriodicPosition::operator +", __FILE__, __LINE__));
}	

return (CxPeriodicPosition(getXPosition() + rhs.getXPosition(),
			   getYPosition() + rhs.getYPosition(),
			   m_xCellSize, m_yCellSize));


}
/** No descriptions */
double CxPeriodicPosition::vabs()const{
	double xTemp = getXPosition();
	double yTemp = getYPosition();
	/*

		Check wether they are on the "right" side of the cell
	*/
	if (xTemp > m_xCellSize*0.5)
	{
		xTemp = xTemp - m_xCellSize;
	}
	if ( yTemp > m_yCellSize*0.5)
	{
		yTemp = yTemp - m_yCellSize;
	}
	double retVal = sqrt ((yTemp*yTemp)+(xTemp * xTemp));
	glLogger.debug("CxPeriodicPosition::vabs returns (%f)", retVal);
	return (retVal);
}
/** Ctor fore conversion of CxPosition to CxPeriodicPosition */
 CxPeriodicPosition::CxPeriodicPosition( const CxPosition &rhs, double aspectRatio){
if (aspectRatio <= 0.0)
{
	throw (CxBadValueError(__FILE__, __LINE__, " Neagtive aspect ratio here!"));
}
m_xCellSize=1.0;
m_yCellSize= aspectRatio;
xPos = mapXPos(rhs.getXPosition());
yPos = mapYpos(rhs.getYPosition()); 
}


/** Subtraktion  operator */
const CxPeriodicPosition CxPeriodicPosition::operator - (const CxPeriodicPosition & rhs){

if ( m_xCellSize != rhs.get_xCellSize() )
{
  glLogger.error(" My xcellSize =(%f), rhs xcellIze =(%f)",m_xCellSize, rhs.get_xCellSize()); 
		throw (CxErrors("veftor not of same cell"));
}
return (CxPeriodicPosition( mapXPos(getXPosition() - rhs.getXPosition()),
		mapYpos(getYPosition() - rhs.getYPosition()), m_xCellSize, m_yCellSize)) ;

}

int CxPeriodicPosition::getNearestPosition(std::vector<CxPeriodicPosition> vPositions) 
{
 glLogger.info("Entering getNearestPosition wiuth positions (%f), (%f)",
		xPos, yPos);
 int retVal = 0;
 double minimalDistance = m_yCellSize; // There must be a start value
    
    if (vPositions.size() == 0)
    {

        return (retVal);
    }
    
    for (unsigned int vIndex = 0; vIndex < vPositions.size(); vIndex++)
    {
      CxPeriodicPosition differenzVec = (*this) - vPositions[vIndex];
      
      double tempDistance =  ((*this)-vPositions[vIndex]).vabs();
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


CxPeriodicPosition & CxPeriodicPosition::operator =(const CxPeriodicPosition &rhs)
{
  if (this != &rhs)
    {
      glLogger.debug("CxPeriodicPosition::operator =.  xCellSize=(%f), (%f)",
		     rhs.get_xCellSize(), rhs.get_yCellSize());
     
      CxPosition::operator=(rhs);
      m_xCellSize = rhs.get_xCellSize();
      m_yCellSize = rhs.get_yCellSize();
    }
  return (*this);

}
