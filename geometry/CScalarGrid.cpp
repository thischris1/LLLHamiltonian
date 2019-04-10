#include <geometry/CScalarGrid.h>
#include <LLLlib/LLLlib.h>
#include <math.h>

#include <iostream>

CScalarGrid::CScalarGrid():
CGrid()
{
}

CScalarGrid::CScalarGrid(float xSize, float ySize , int xMeshCount, int yMeshCount):
CGrid(xSize, ySize, xMeshCount, yMeshCount)
{
	
}

CScalarGrid::~CScalarGrid()
{
}

double CScalarGrid::getValueAtIndex(const int &xIndex, const int & yIndex)const
{	
	double xAufpunkt = 0.0;
	double yAufpunkt = 0.0; 
	if (!CGrid::meshIndexToCoordinate(xIndex, yIndex,xAufpunkt, yAufpunkt))
			{
				return (0);
			}
			return( getValue(CxPosition(xAufpunkt, yAufpunkt)));	
}


double CScalarGrid::getValueAtPosition(const CxPosition &n_pos)const
{
	return (getValue(n_pos));;	
	
}
/*! calculate a correlation
 */
CDistanceDistribution    CScalarGrid::getCorrelation(const int steps, double eps)
{
	float length = sqrt(CGrid::getXLength()* CGrid::getXLength() + CGrid::getYLength()* CGrid::getYLength());
	
	if (steps <10) 
	{
		CDistanceDistribution retVal(0,0.0f);
		return (retVal);
	}
	// calculate needed precision, should be 1/100 of the maximal spread
	double 
		max = 0.0, 
		min = 0.0;
	//int complete =  CGrid::getxSteps() * CGrid::getySteps();
		if (!getExtrema(max,min))
		{
			ERRORTHROW("error in getExtrema");
		}
	/*
		eps = (fabs(max-min)/100.0);
		if (eps < 1e-06)
		{
			eps = 1e-06;
		}
		*/
	CDistanceDistribution retVal(steps, length); 
	// loop over all points
	int progressCount = 0;
	for (int xCor = 0; xCor < CGrid::getxSteps(); ++xCor) {
		for (int yCor = 0; yCor < CGrid::getySteps(); ++yCor) {
			double targetValue = getValueAtIndex(xCor, yCor);
			 /*
			  * Inner loop over all points to find the same value again
			  */
			  CxPosition aufPunkt=indexToPosition(xCor, yCor);
			  glLogger.info("Aufpunkt = (%f),(%f)", aufPunkt.getXPosition(), aufPunkt.getYPosition());
			for (int iXpos = 0; iXpos < CGrid::getxSteps(); ++iXpos) {
				for (int iYPos = 0; iYPos < CGrid::getySteps(); ++iYPos) {
					double tempVal = getValueAtIndex(iXpos, iYPos);
					if (fabs(tempVal - targetValue) < eps)
					{
						/*
						 * FOUND
						 * calculate distance to Aufpunkt
						 * increment vector at point
						 */	
						 CxPosition found =indexToPosition(iXpos,iYPos);
						glLogger.info("gefunden bei  = (%f),(%f)", found.getXPosition(), found.getYPosition());
						 float distance = (found - aufPunkt).vabs();
						 retVal.addValue(distance);
						 glLogger.info("Found the same value at a distance (%f)", distance);
						 
					}
					
					
				}
				progressCount++;
			}
			
		}
		glLogger.error("Line (%i) of (%i)", xCor,  CGrid::getxSteps());
	}
	// save the distribution temporarily
	glLogger.error("Writing distribution");
	std::ofstream out("./testDistribution.dat");
	retVal.writeToStream(out);
	out.close();
	 	return (retVal);
}

CxPosition CScalarGrid::findMaximumPosition(void)const
{
	double retVal = -1e36;
	int xPos = 0;
	int yPos = 0;
	for (int xIndex = 0; xIndex < getxSteps(); ++xIndex) {
		for (int yIndex = 0; yIndex < getySteps(); ++yIndex) {
			double tempVal = getValue(indexToPosition(xIndex, yIndex));
			if (retVal < tempVal)
			{
				// new minimum
				retVal = tempVal;
				xPos = xIndex;
				yPos = yIndex;
			}			
		}
		
	}
	
	glLogger.debug("Maximum (%f) at x= %f, y = %f",retVal, xPos, yPos);
	return(indexToPosition(xPos,yPos));

	
}

CxPosition CScalarGrid::findMinimumPosition(void)const
{
	double retVal = 1e36;
	int xPos = 0;
	int yPos = 0;
	for (int xIndex = 0; xIndex < getxSteps(); ++xIndex) {
		for (int yIndex = 0; yIndex < getySteps(); ++yIndex) {
			if (retVal > getValue(indexToPosition(xIndex, yIndex)))
			{
				// new minimum
				retVal = getValue(indexToPosition(xIndex, yIndex));
				xPos = xIndex;
				yPos = yIndex;
			}			
		}
		
	}
	glLogger.debug("Minimum %f at x= %f, y = %f",retVal, xPos, yPos);
	 
	return (indexToPosition(xPos,yPos));
}

bool CScalarGrid::dumpContent(void) const
{
	int m_size = getxSteps()*getySteps();
	std::cout << " Total Size = "<< m_size<<" xSize = "<< getxSteps() << " ySize = " << getySteps()<<std::endl;
	for (int rowIndex = 0; rowIndex < getxSteps(); ++rowIndex) {
		for (int colIndex = 0; colIndex < getySteps(); ++colIndex) {
			CxPosition tempPos = indexToPosition(rowIndex, colIndex);
			std::cout << "("<< rowIndex<<"), (" << colIndex << ") = (" << getValue(tempPos) <<std::endl;
			
		}
		
	}
	return (true);
}

double CScalarGrid::getMaximum(void)const
{
	
	double 
		max = 0.0, 
		min = 0.0;
	if (getExtrema(max,min))
	{
		return (max);
	}
	else
	{
		//ERROR
		return (-1e36);
	}
}
double CScalarGrid::getMinimum(void)const
{
	
	double 
		max = 0.0, 
		min = 0.0;
	if (getExtrema(max,min))
	{
		return (min);
	}
	else
	{
		//ERROR
		return (1e36);
	}
}

bool CScalarGrid::fillGrid(IValueProvider *m_provider)
{
	for (int xIndex = 0; xIndex < getxSteps(); ++xIndex) {
		for (int yIndex = 0; yIndex < getySteps(); ++yIndex) {
			CxPosition tempPos = indexToPosition(xIndex,yIndex);
			double tempVal = m_provider->getValue(tempPos);
			
			setValue(tempPos,tempVal);
			
		}
		
	}	
	return (true);
}


