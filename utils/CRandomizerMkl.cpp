#include "CRandomizerMkl.h"
//#include <mkl/mkl_vsl.h>
#include <time.h>
CRandomizerMkl * CRandomizerMkl::m_instance = 0;

CRandomizerMkl::CRandomizerMkl()
{
	
	vslNewStream(&stream,VSL_BRNG_MCG31,2);

}

CRandomizerMkl::~CRandomizerMkl()
{
  vslDeleteStream(&stream);
  
}
int CRandomizerMkl::iRand(const int start, const int end)
{
	int retVal = 0;
	if (start == end)
	{
		
	}
	else 
	{
	
		retVal = retVal + start;
	}
	return(retVal);	
	
}

float CRandomizerMkl::fRand(const float start, const float end)
{
	
	
	return ((float) dRand(start,end));	
}

double CRandomizerMkl::dRand(const double start, const double end)
{
	double retVal = 0.0;
	
	retVal = retVal * (end - start) + start;
	return (retVal);	
}

CRandomizerMkl * CRandomizerMkl::getInstance()
{
	if (m_instance == 0)
    {
      m_instance = new CRandomizerMkl();
    }
  return m_instance;
}

void CRandomizerMkl::randomize(int seed)
{
	if (seed == 0)
	{
		time_t secs  = time(NULL);


	}
	else
	{

	}


}
