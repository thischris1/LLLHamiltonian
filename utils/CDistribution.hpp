#ifndef CDISTRIBUTION_H
#define CDISTRIBUTION_H

#include <vector>
#include <fstream>

class CDistribution {
 public:
  CDistribution();
  CDistribution (int);
  CDistribution (const CDistribution & rhs);
  CDistribution & operator =(const CDistribution &rhs);
  virtual ~CDistribution();
  virtual bool addValue(const float);
  virtual bool addValue (const double);
  virtual bool addValueAtIndex(const int index);
  virtual bool writeToStream(std::ostream &) const;
  int getSize(void) const;
  std::vector <int> getHistogram(void)const;
  
 protected:
  bool init (void);
  virtual float calculateArea(float & distance) const;

 private:

 
  int N_STEP;

 std::vector<int> m_Func;







};




#endif


