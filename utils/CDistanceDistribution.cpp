#include <utils/CDistanceDistribution.h>
#include <math.h>
#include <LLLlib/LLLlib.h>

CDistanceDistribution::CDistanceDistribution(void)
: CDistribution(100), xMax(1.0f)
,xMin(0.0f)
{

}

CDistanceDistribution::~CDistanceDistribution(void)
{
}

CDistanceDistribution::CDistanceDistribution(int gridSize , float  maxDistance, float minDistance)
: CDistribution(gridSize),
xMax(maxDistance),
xMin(minDistance)
{
	if (maxDistance < 0.0)
	{
		
		ERRORTHROW("Wrong order of min and max");
	}
} 

CDistanceDistribution::CDistanceDistribution(int gridSize , float  maxDistance):
CDistribution(gridSize),
xMax(maxDistance),
xMin(0.0f)
{
	if (maxDistance < 0.0f)
	{
		ERRORTHROW("Wrong order of min and max");
	}
	
}


bool CDistanceDistribution::addValue(const float fIndex)
{
  glLogger.info("Distance Distro: Add distance (%f)", fIndex);
  return(CDistribution::addValueAtIndex(calcIndex(fIndex)));
}

int CDistanceDistribution::calcIndex(const float fIndex )
{
	if (fIndex <= xMax && fIndex >= xMin)
	{
		float stepSize = (xMax-xMin)/getSize();

		int tempVal = (int)floor(fIndex/stepSize);
		return (tempVal);
	}
	return (-1);
}

bool CDistanceDistribution::writeToStream(std::ostream & outStream)
{
  if (!outStream.good())
    {
      return (false);
    }
  float stepSize = (xMax -xMin) /getSize();
  float xVal = xMin;
  std::vector <int> result = CDistribution::getHistogram();
  for (int index = 0; index < getSize(); index ++) 
    {
      outStream << xVal<<" "<< result[index]<<std::endl;
      xVal = xVal+ stepSize;
    }
  return (true);
}

CDistanceDistribution::CDistanceDistribution(const CDistanceDistribution & rhs): CDistribution(rhs)
{
	if (this != &rhs)
	{
		xMax = rhs.getMax();
		xMin = rhs.getMin();

	}
}




float CDistanceDistribution::getMax(void) const
{
	return (xMax);
}

float CDistanceDistribution::getMin(void) const
{
	return (xMin);
}

float CDistanceDistribution::calculateArea(float & distance) const
{
  if (distance <= 0)
    {
      return (0.0f);
    }
  float stepSize = (xMax -xMin) /getSize();
  float retVal = stepSize*stepSize + 2* stepSize*distance;
  retVal = M_PI*retVal*0.25f;
  if (xMax == xMin)
    {
      // square case

      if (distance < xMax)
	{
	 
	  return (retVal);
	}
    }
  else 
    {
      // calculate opening angle 
      
      float tempAngle = acos(distance/xMax);
      tempAngle = 1 - 4*tempAngle/M_PI;
      retVal = tempAngle * retVal;
    }
  return (retVal);
}
CDistanceDistribution & CDistanceDistribution::operator =(const CDistanceDistribution & rhs)
{
	if (&rhs != this)
	{
		CDistribution::operator =(rhs);
		xMax = rhs.getMax();
		xMin = rhs.getMin();	
		
	}
	return (*this);
	
}

