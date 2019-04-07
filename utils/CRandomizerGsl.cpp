#include "CRandomizerGsl.h"
#include <gsl/gsl_rng.h>
#include <time.h>
CRandomizerGsl * CRandomizerGsl::m_instance = 0;

CRandomizerGsl::CRandomizerGsl()
{
	
	
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
}

CRandomizerGsl::~CRandomizerGsl()
{
	gsl_rng_free(r);
	r = 0;
}
int CRandomizerGsl::iRand(const int start, const int end)
{
	int retVal = 0;
	if (start == end)
	{
		retVal = gsl_rng_get(r);
	}
	else 
	{
		retVal = gsl_rng_uniform_int(r,end-start);
		retVal = retVal + start;
	}
	return(retVal);	
	
}

float CRandomizerGsl::fRand(const float start, const float end)
{
	
	
	return ((float) dRand(start,end));	
}

double CRandomizerGsl::dRand(const double start, const double end)
{
	double retVal = 0.0;
	retVal = gsl_rng_uniform(r);
	retVal = retVal * (end - start) + start;
	return (retVal);	
}

CRandomizerGsl * CRandomizerGsl::getInstance()
{
	if (m_instance == 0)
    {
      m_instance = new CRandomizerGsl();
    }
  return m_instance;
}

void CRandomizerGsl::randomize(int seed)
{
	if (seed == 0)
	{
		time_t secs  = time(NULL);
		gsl_rng_set(r,secs);

	}
	else
	{
		gsl_rng_set(r,seed);
	}


}
