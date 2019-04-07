#ifndef RANDOMIZER_H
#define RANDOMIZER_H
/*
 APPLICATION INCLUDES
*/
#include <geometry/CxPosition.hpp>
#include <geometry/CxPeriodicPosition.hpp>
// SYSTEM INCLUDES

#include <vector>

/*!

\file has the declaration of the randomizer class

*/

extern "C" {
  float ran0(long &);
  float ran1(long &);
}


class CRandomizer {

 public:
  static CRandomizer * getInstance(void);

  virtual int iRand(const int start = 0, const int  end = 0 );
  virtual float fRand(const float start = 0.0f, const float  end = 1.0f);
  virtual double dRand(const double start =0.0 , const double  end = 1.0);

  int randomPositions(std::vector<CxPosition> & positions);
  int randomPositions(std::vector<CxPeriodicPosition> & positions);
   virtual void randomize(int seed = 0);
   CRandomizer();  
  virtual ~CRandomizer();
 protected:


 private:
    long idum;

  CRandomizer (CRandomizer &);
  static CRandomizer * instance;

};




#endif 
