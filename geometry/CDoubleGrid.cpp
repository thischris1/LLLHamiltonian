#include "CDoubleGrid.h"
#include <list>
#include <iostream>
CDoubleGrid::CDoubleGrid(float xSize, float ySize , int xMeshCount, int yMeshCount)
: CScalarGrid(xSize, ySize, xMeshCount, yMeshCount)
{
	initialize();	
	
}


CDoubleGrid::CDoubleGrid()
{
}

CDoubleGrid::~CDoubleGrid()
{
	m_values.clear();
}

double CDoubleGrid::getValue(const CxPosition &n_pos)const
{
	int index = positionToIndex(n_pos);
	return (m_values[index]);
	
}
bool CDoubleGrid::initialize(double newVal)
{
	m_values = std::vector<double> (getxSteps()*getySteps(), newVal);
	return (true);
}

bool CDoubleGrid::setValue (const CxPosition &n_pos, double newValue)
{
	int index = positionToIndex(n_pos);
	m_values[index] = newValue;
	return (true);	
	
}	
	
bool CDoubleGrid::getExtrema(double &max, double & min) const
{
	double tempMax= -1e36;
	double tempMin = 1e36;
	std::vector<double>::const_iterator it;
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
	
	
bool CDoubleGrid::dumpFullContent(void)
{
	std::vector<double>::const_iterator it;
	it  = m_values.begin();
	int count = 0;
	while (it != m_values.end())
	{
		std::cout << "Index ("<<count<<") = " << *it <<std::endl;
		count++;
		it++;
	} 
	
return (true);	
}
	
