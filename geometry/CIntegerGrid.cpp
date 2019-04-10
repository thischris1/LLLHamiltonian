#include "CIntegerGrid.h"
#include <math.h>
#include <utils/logger.hpp>

CIntegerGrid::CIntegerGrid()
{
}

CIntegerGrid::~CIntegerGrid()
{
	/*
	 * Loop over all rows in matirx, delete them
	 */
	 m_values.clear();
}

CIntegerGrid::CIntegerGrid(float xSize, float ySize , int xMeshCount, int yMeshCount)
:CScalarGrid(xSize, ySize, xMeshCount, yMeshCount)
{
	initialize();	
	
}
bool CIntegerGrid::initialize()
{
	/*
	 * Allocate array here
	 */	
	 
	 int totalSize =getxSteps() * getySteps();
	 glLogger.debug("Creating a vector of size %i", totalSize); 
	m_values = std::vector<int>(totalSize,0);
	 
	 return (true);
}

double CIntegerGrid::getValue(const CxPosition &n_pos) const
{
	// convert to mesh point
	int index = positionToIndex(n_pos);
	return (m_values[index]);
	
	
}

bool CIntegerGrid::setValue(const CxPosition &n_pos, double newValue)
{
	int index = positionToIndex(n_pos);
	m_values[index] = static_cast<int>(round(newValue));
	return (true); 	
	
}

bool CIntegerGrid::getExtrema(double &max, double & min) const
{
	double tempMax= -1e36;
	double tempMin = 1e36;
	std::vector<int>::const_iterator it;
	it = m_values.begin();

	while (it != m_values.end())
	{
		if (*it < tempMin)
		{
			tempMin = *it;
		}
		if (*it > tempMax)
		{
			tempMax = *it;	
		}
		
		
		it++;
	} 
	max = tempMax;
	min = tempMin;	

	return (true);	
}


bool CIntegerGrid::dumpFullContent(void)
{
	std::vector<int>::iterator it;
	it = m_values.begin();
	int count = 0;
	while (it != m_values.end())
	{
		std::cout << "Index ("<<count<<") = " << *it <<std::endl;
		count++;
		it++;
	} 
	
return (true);	
}

