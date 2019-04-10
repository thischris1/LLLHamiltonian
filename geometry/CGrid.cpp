
#include "CGrid.h"
#include <math.h>
#include <utils/logger.hpp>
CGrid::CGrid():
m_xLength(1.0),
m_yLength(1.0),
m_xSteps(100),
m_ySteps(100),
m_xMeshSize(0.01),
m_yMeshSize(0.01)
{

}

CGrid::CGrid(double xSize, double ySize, int xSteps, int ySteps):
m_xLength(xSize),
m_yLength(ySize),
m_xSteps(xSteps),
m_ySteps(ySteps),
m_xMeshSize(0.01),
m_yMeshSize(0.01)
{
	if ( (m_xLength < 0) || (m_yLength < 0 ) || (m_xSteps <= 1) || (m_ySteps <= 1))
	{
		// throw something here	
	}
	m_xMeshSize = m_xLength / m_xSteps;
	m_yMeshSize = m_yLength / m_ySteps; 
}

CGrid::~CGrid()
{
	// reset values
}
/*
 * 
 */
int CGrid::getNearestXIndex(const double & xPoint) const {

	if ( (xPoint < 0.0) || (xPoint > m_xLength) )
	{
		return (-1);	
	} 	
	int retVal = 0;
	retVal = static_cast<int>(round(xPoint/m_xMeshSize));
	return retVal;
	
}

int CGrid::getNearestYIndex(const double & yPoint) const {

	if ( (yPoint < 0.0) || (yPoint > m_yLength) )
	{
		return (-1);	
	} 	
	int retVal = 0;
	retVal = static_cast<int>(round(yPoint/m_yMeshSize));
	return retVal;
	
}

double CGrid::getNearestXPoint(const double & xPoint) const {
	return (getNearestXIndex(xPoint)*m_xMeshSize);
}

double CGrid::getNearestYPoint(const double & yPoint) const {
	return (getNearestYIndex(yPoint)*m_yMeshSize);
}

int CGrid::mapIndex (int xPoint, int yPoint) const
{
	if ( (xPoint < 0) || (yPoint < 0) || (xPoint > m_xSteps) || (yPoint > m_ySteps) )
	{
		return (-1);
	}
	int retVal = yPoint*m_xSteps + xPoint;
	glLogger.debug("mapIndex transforms x = (%i) and y = (%i) to (%i)",xPoint, yPoint,retVal); 
	return (retVal);
}
	
int CGrid::mapDoubleToIndex (double xPoint, double yPoint) const
{
	int xTemp = getNearestXIndex(xPoint);
	if (xTemp < 0 )
	{
		return -1;
	}
	int yTemp = getNearestYIndex(yPoint);
	if (yTemp < 0 )
	{
		return -1;
	}	
	return (mapIndex(xTemp, yTemp));
}

int CGrid::positionToIndex(const CxPosition &nPos)const
{
	
	return (mapDoubleToIndex(nPos.getXPosition(), nPos.getYPosition()));
}
int CGrid::meshIndexToCoordinate(int x, int y, double & n_xPos, double & n_yPos)const 
{
	if ( (x < 0) || (y < 0 ) || (x > m_xSteps) || (y > m_ySteps))
	{ 
		n_xPos = -1.0;
		n_yPos = -1.0;
		return (0);
	}	
	n_xPos = x*m_xMeshSize;
	n_yPos = y*m_yMeshSize;
	return (1);
}


CxPosition CGrid::indexToPosition(int rowIndex, int colIndex) const
{
double xTemp, yTemp;
if (meshIndexToCoordinate(rowIndex, colIndex,xTemp,yTemp))
{
	
	return (CxPosition(xTemp, yTemp));
}
else 
{
	// ERROR HERE
	glLogger.error("Cant map");
	return (CxPosition(-1,-1));	
}
	
}


