#ifndef CDISTANCEDISTRIBUTION_H
#define CDISTANCEDISTRIBUTION_H
#include <utils/CDistribution.hpp>
#include <fstream>
class CDistanceDistribution :
	public CDistribution
{
public:
	CDistanceDistribution(void);
	~CDistanceDistribution(void);
private:
	// The maximal value 
	float xMax;
	float xMin;
public:
	CDistanceDistribution(int gridSize , 
			      float  maxDistance, 
			      float minDistance);
	CDistanceDistribution(int gridSize , float  maxDistance);
	bool addValue(const float fIndex);
	int calcIndex(const float fIndex);
	bool writeToStream(std::ostream & outStream);
	CDistanceDistribution(const CDistanceDistribution & rhs);
	CDistanceDistribution & operator =(const CDistanceDistribution & rhs);
	float getMax(void) const;
	float getMin(void) const;
	float calculateArea(float & distance) const;
};
#endif
