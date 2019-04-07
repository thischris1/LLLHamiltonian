#include <utils/CRandomizer.h>
#include <utils/logger.hpp>
#include <time.h>


CRandomizer::CRandomizer()
{
  idum = 34156283;

}

CRandomizer::~CRandomizer()
{
	/*
	 * Nothing so far
	 *  
	 */
	 instance = 0;
}

CRandomizer::CRandomizer(CRandomizer &)
{


}

CRandomizer * CRandomizer::getInstance()
{
  if (instance == 0)
    {
      instance = new CRandomizer();
    }
  return instance;
}


int CRandomizer::iRand(const int  start, const int end)
{

  return 1;

}

float CRandomizer::fRand(const float  start, const float end)
{

  float tempVal = ran1(idum);
  if (0 == start && end == 1)
    {
      return (tempVal);
    }
  else 
    {
      float differenz = end -start;
      tempVal = start + differenz * ran1(idum);
      return (tempVal);
    }

}



double CRandomizer::dRand(const double   start, const double end) 
{

  return ((double)fRand((float)start, (float) end));


}

int CRandomizer::randomPositions(std::vector<CxPosition> &positions)
{
  int oldSize = (int) positions.size();
  for (int index=0; index <oldSize; index++)
    {
      positions[index].setXPosition(fRand());
      positions[index].setYPosition(fRand());
      CxPosition tempPosition(fRand(), fRand());\
     glLogger.info ("idum is (%d)", idum);
      positions[index]= tempPosition;
      glLogger.info("Vector holds now (%f), (%f) at index (%d)", 
		    positions[index].getXPosition(), 
		    positions[index].getYPosition(), index);
    }
  return (positions.size());
}
int CRandomizer::randomPositions(std::vector<CxPeriodicPosition> &positions)
{
  int oldSize = (int) positions.size();
  for (int index=0; index <oldSize; index++)
    {
      positions[index].setXPosition(fRand());
      positions[index].setYPosition(fRand());
      CxPeriodicPosition tempPosition(fRand(), fRand());
     glLogger.info ("idum is (%d)", idum);
     positions[index]= CxPeriodicPosition(tempPosition);
      glLogger.info("Vector holds now (%f), (%f) at index (%d)", 
		    positions[index].getXPosition(), 
		    positions[index].getYPosition(), index);
    }
  return (positions.size());
}



#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long &idum)
{
	long k;
	float ans;

	idum ^= MASK;
	k=(idum)/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	ans=AM*(idum);
	idum ^= MASK;
	return ans;
}



#undef MASK



#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long & idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (idum <= 0 || !iy) {
		if (-(idum) < 1) idum=1;
		else idum = -(idum);
		 for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=(idum)/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


void CRandomizer::randomize(int seed)
{
	if (seed == 0)
	{
		// make seed the seconds and minutes	
		time_t ticks = time(0);
		tm * m_time = gmtime(&ticks);
		int count = m_time->tm_sec + m_time->tm_mday + m_time->tm_hour;
		for (int index = 0; index < count; ++index) {
			float temp = fRand(0.0f, 1.0f);
			temp = temp+1.0;
			
		}
	}	
	
}

 CRandomizer * CRandomizer::instance = 0;
