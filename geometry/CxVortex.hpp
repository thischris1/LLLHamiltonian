#ifndef CXVORTEX_HPP
#define CXVORTEX_HPP

#include <geometry/CxPeriodicPosition.hpp>

class CxVortex : public CxPeriodicPosition
{
public:
  CxVortex ();

  CxVortex(double xPos, double yPos, int windingNumber, 
	   double xSize = 1.0, double ySize = 1.0);

  CxVortex (CxPeriodicPosition, int);
  int getWindingNumber(void) const;
  
  void setWindingNumber (int);
  /**Assignement opeator */
  CxVortex operator =(const CxVortex & );
  /** No descriptions */
   CxVortex (const CxVortex & rhs);
  const CxVortex  operator-(const CxVortex &rhs);
private:
  int windNumber; 



};











#endif
