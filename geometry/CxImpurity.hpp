
/*!
	\file CxImpurity.hpp
        \brief Declaration of class CxImpurity
	\author  Christian Mueller
	\date    07 May 03

	
*/ 
#ifndef CXIMPURITY_HPP
#define CXIMPURITY_HPP
#include <geometry/CxPosition.hpp>
/*!
  \class CxImpurity
  \brief encapsulates an gaussian shaped impurity in a a x b cell with strength, sigmax and sigma y 
*/

class CxImpurity : public CxPosition
{
public:
    //! Default ctor
  CxImpurity();
    //! Standard Constructor
  CxImpurity( double xPos, double yPos, double Nsigmax, double Nsigmay, double Nstrength);
 
 CxImpurity (CxPosition &pos, double Nsigmax, double Nsigmay, double Nstrength);
  CxImpurity (const CxImpurity &rhs);
/** assignement */
  CxImpurity operator = (const CxImpurity &rhs);
  // 
  virtual ~CxImpurity();
  //! Accessor function for strength
  double getStrength(void) const{ return (strength);};
  //! Accessor method for sigma x,
  double getSigmaX(void) const { return(sigmax);};
  //! Accessor method for sigma y
  double getSigmaY(void) const { return(sigmay);};
  //! Set methods for x,y, sigma
  void setSigmaX(double newSigmax);
  void setSigmaY(double newSigmay);
  void setStrength (double newStrength);
  bool operator == (const CxImpurity &rhs);
  float getPotential(const CxPosition & aufPunkt) const; 

private:
  
  /** Copy ctor */
  double sigmax;
  double sigmay;
  double strength;

};

#endif
