#ifndef C2DDISTRIBUTION_H
#define C2DDISTRIBUTION_H
#include <vector>
#include <fstream>
#include <vector>
#include <geometry/C2DVector.hpp>
#include <geometry/CxPosition.hpp>
#include <geometry/CxPeriodicPosition.hpp>

class C2dDistribution {
public:
  //! default ctor
  C2dDistribution();
  //! Ctor with given size
  C2dDistribution (const int & x_size, const int ySize);
  //! dtor
  virtual ~C2dDistribution();
  //! copy ctor
  C2dDistribution(const C2dDistribution &rhs);
  //! basic addValue method
  bool addValue (const CxPosition & newPos);
  //! add a periodic position here
  bool addValue (const CxPeriodicPosition &newPos);


  //! Write to a stream 
  virtual bool writeToStream(std::ostream & outStream) const;
  //! gets the array vector
  C2DVector<int>  getValues(void) const;

  //! access to size
  unsigned int getXSteps() const {return (N_STEPS_X);}
  unsigned int getYSteps() const {return (N_STEPS_Y);}
private:
  unsigned int N_STEPS_X;
  unsigned int N_STEPS_Y;
  C2DVector <int>  m_Values;
};















#endif
