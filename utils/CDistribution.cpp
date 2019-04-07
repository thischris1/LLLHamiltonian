#include <utils/CDistribution.hpp>
#include <math.h>
#include <LLLlib/LLLlib.h>
CDistribution::CDistribution()
  :N_STEP(100),
   m_Func(std::vector<int>(100,0))
{
  N_STEP =100;
  init();
}



CDistribution::CDistribution(int n_Size)
  :N_STEP(n_Size),
  m_Func(std::vector<int>(n_Size,0))
{


}

CDistribution::CDistribution( const CDistribution & rhs):
  N_STEP(rhs.getSize())
{
	if  (&rhs != this)
	{
  		m_Func=(std::vector<int>(rhs.getHistogram()));
	}

}

CDistribution & CDistribution::operator=(const CDistribution &rhs)
{
  if (&rhs != this)
    {
      N_STEP=rhs.getSize();
      m_Func= rhs.getHistogram();
    }
  return (*this);

}









bool CDistribution::addValue(const float newValue)
{
  return (addValue((double)newValue));
}



bool CDistribution::addValueAtIndex(const int index)
{
	glLogger.info("addValueAtIndex %i, where the initial value is %i", index, m_Func.at(index));
	try {
  m_Func.at(index) =  m_Func.at(index)+1;
  }
  catch(const exception& e)
	{
		ERRORTHROW( e.what());
	}
//  m_Func[tempIndex] =  m_Func[tempIndex] + 1;
	glLogger.info("Value after insertion");
  return (true);
}


bool CDistribution::addValue(const double newValue)
{
  glLogger.info("CDistribution adds (%f)", newValue);
  int tempIndex = (int) floor(newValue*N_STEP);
 
  try {
  m_Func.at(tempIndex) =  m_Func.at(tempIndex) + 1;
  }
  catch(const exception& e)
	{
		ERRORTHROW( e.what());
	}
//  m_Func[tempIndex] =  m_Func[tempIndex] + 1;
  return (true);

}

bool CDistribution::init(void)
{
  m_Func= std::vector<int>(N_STEP,0);
  return (true);
}


CDistribution::~CDistribution()
{
  m_Func.clear();
}

bool CDistribution::writeToStream(std::ostream &outStream) const
{
 
  
  for (int index = 0; index < N_STEP; index++)
    {
      float tempVal = (float) index;
      outStream << index << " \t "<< m_Func[index] << " \t " << calculateArea(tempVal) << std::endl;
    }
  return (true);
}

int  CDistribution::getSize(void) const 
{
  return (N_STEP);
}


std::vector<int> CDistribution::getHistogram(void)const
{

  return (m_Func);
}

float CDistribution::calculateArea(float & distance) const
{
  float xMax = sqrt(0.5*N_STEP*N_STEP); 
 glLogger.debug(" CDistribution::calculateArea (%f)", distance); 
 if (distance < 0)
    {
      
      return (0.0f);
    }
  float stepSize = (N_STEP) /getSize();
  float kreis1 = M_PI*0.25*distance*distance;
  float kreis2 =  M_PI*0.25*(distance+stepSize)*(distance+stepSize);
  float retVal = kreis2-kreis1;

  glLogger.debug("retVal (kreisring flaeche %f", retVal);
 // assuming  square geometry

      if (distance < xMax)
	{
	 
	  return (retVal);
	}
      else 
	{
	  // calculate opening angle 
	  glLogger.debug("Angle calculation %d", distance);
	  float tempAngle = acos(xMax/distance);
	  tempAngle = 1 - 4*tempAngle/M_PI;
	  retVal = tempAngle * retVal;
	}
  return (retVal);
}
