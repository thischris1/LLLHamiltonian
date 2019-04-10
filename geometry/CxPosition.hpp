/*!
  \file CxPosition.hpp
  \brief Contains declaration of class CxPosition
  \author Christian Mueller
  \date 07 May 03
*/

#ifndef CXPOSITION_HPP
#define CXPOSITION_HPP
#include <vector>

/*!
  \class CxPosition
  \brief encapsulates the 2-d geometry
 */

class CxPosition
{
public:
   //! Standard ctor
  CxPosition(); 
    //! Constructor 
  CxPosition(double newX , double newY);   
    //! getter method
  double getXPosition(void) const ;
    //! getter method
  double getYPosition(void) const ;
    //! Setter method
  void setXPosition (double newXpos) ;
  void setYPosition (double newYpos) ;
  //! Compares two positions

  bool compare (const CxPosition &)const;
    //! Setter for epsilon
  void setEpsilon(double);
    //! getter for epsilon
  double getEpsilon(void) const;
  //! Compares two Objects (see compare())
  virtual bool operator == (const CxPosition &rhs);
  //! Vector subtraction
  
   const CxPosition   operator - (const CxPosition &);

   const CxPosition   operator + (const CxPosition &);
  
  CxPosition &operator = (const CxPosition & rhs);

  void multiplyVector (double );
//! Absolute value  
  virtual double vabs(void) const;

  virtual double vabs (CxPosition) const;
  
    //! Get the nearest position from the vPositions
   virtual int  getNearestPosition(std::vector <CxPosition> vPositions);

  /** No descriptions */
  virtual  ~CxPosition();
  /** Copy ctor */
   CxPosition(const CxPosition &rhs);


protected:
  double xPos;
  double yPos;
  double epsilon; //! Allows for "fuzzy" comparison of positions
  
 private:
 

};

/*!
  \class class CxBoundPosition: public CxPosition
  \brief encapsulates the FQH geometry i.e. a cell with length 1 and width 1*aspectratio. Positions outside this geometry will be translated back into this unit cell 
  

*/
class CxBoundPosition: public CxPosition
{
public:
  // Standrad constructor
  CxBoundPosition();
  //! Destructor
  virtual ~CxBoundPosition();

  //! Constructor to be used 
  CxBoundPosition (double xPos, double yPos, double aspectRatio);
  //! Setter function for aspectratio
  void setAspectRatio(double);
protected:
  double aspectRatio;
};










#endif 















