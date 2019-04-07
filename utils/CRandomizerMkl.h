#ifndef CRANDOMIZERMKL_H_
#define CRANDOMIZERMKL_H_

#include "CRandomizer.h"
#include <mkl.h>

class CRandomizerMkl : public CRandomizer
{
public:
	CRandomizerMkl();
	virtual ~CRandomizerMkl();
	static CRandomizerMkl * getInstance(void);

  int iRand(const int start = 0, const int  end = 0 );
  float fRand(const float start = 0.0f, const float  end = 1.0f);
  void randomize(int seed = 0);
  double dRand(const double start =0.0 , const double  end = 1.0);
	private:
  VSLStreamStatePtr stream;
  static CRandomizerMkl *m_instance;


	
};

#endif /*CRANDOMIZERMKL_H_*/
