#ifndef CRANDOMIZERGSL_H_
#define CRANDOMIZERGSL_H_

#include "CRandomizer.h"
#include <gsl/gsl_rng.h>

class CRandomizerGsl : public CRandomizer
{
public:
	CRandomizerGsl();
	virtual ~CRandomizerGsl();
	static CRandomizerGsl * getInstance(void);

  int iRand(const int start = 0, const int  end = 0 );
  float fRand(const float start = 0.0f, const float  end = 1.0f);
  void randomize(int seed = 0);
  double dRand(const double start =0.0 , const double  end = 1.0);
	private:
	
	static CRandomizerGsl *m_instance;
	const gsl_rng_type *T;
	gsl_rng *r;
	
};

#endif /*CRANDOMIZERGSL_H_*/
