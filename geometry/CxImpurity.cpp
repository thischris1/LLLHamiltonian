/*!
  \file CxImpurity.cpp
  \brief Implementation of class CxImpurity
  \author Christian Mueller

*/
#include <geometry/CxImpurity.hpp>
#include <utils/logger.hpp>

CxImpurity::CxImpurity(): CxPosition(0,0)
 
{
  sigmax = 0;
  sigmay = 0;
  strength = 0;

}

CxImpurity::CxImpurity (CxPosition &pos, double Nsigmax, double Nsigmay, double Nstrength):
CxPosition(pos),
sigmax(Nsigmax),
sigmay(Nsigmay),
strength(Nstrength)
{
	
}


CxImpurity::CxImpurity(double xPos, double yPos,
		       double Nsigmax, double Nsigmay, 
		       double Nstrength)
  :CxPosition(xPos, yPos)
{

  sigmax =  Nsigmax;
  sigmay =  Nsigmay;
  strength = Nstrength; 


}
/*!
  \fn void CxImpurity::setSigmaX(double newSigmax)
  \brief setter for \sigma_x
  \param  newSigmax the new x-shape of the gaussian
*/

void CxImpurity::setSigmaX(double newSigmax)
{
  if (newSigmax < 0 )
    {
      sigmax= 0;
      return;
      //! \todo throw an error here
    }
  sigmax = newSigmax;
}


/*!
 \fn void CxImpurity::setSigmaY(double newSigmay)
 \brief sets the y-shape of the gaussian
 \param  newSigmay the new Sigma
*/

void CxImpurity::setSigmaY(double newSigmay)
{
  if (newSigmay <0)
    {
      sigmay = 0;
      return;
    }
  sigmay = newSigmay;

}

/*!
  \fn void CxImpurity::setStrength (double newStrength)
  \brief setter for strength
  \param newStrength the new strength of the impurity
*/

void CxImpurity::setStrength (double newStrength)
{
  
  strength = newStrength;


}
/** Copy ctor */
 CxImpurity::CxImpurity(const CxImpurity & rhs):
 CxPosition(rhs)
 {
 	sigmax	= rhs.getSigmaX();
	sigmay = rhs.getSigmaY();
	strength = rhs.getStrength();
}
/** No descriptions */
CxImpurity CxImpurity::operator = (const CxImpurity &rhs){
if ( &rhs != this){
	CxPosition::operator =(rhs);
	sigmax	= rhs.getSigmaX();
	sigmay = rhs.getSigmaY();
	strength = rhs.getStrength();
	}
	return (*this);
}

bool CxImpurity::operator ==(const CxImpurity &rhs)
{
	if ( (fabs(rhs.getSigmaX() - sigmax) > 1e-05))
	{
		glLogger.debug("Different positions");
		return (false);
	}
	else if ( fabs(rhs.getSigmaY() -sigmay) > 1e-05 ) 
	{
		glLogger.debug("Different positions");
		return (false);	
	}
	else if ( fabs(rhs.getStrength() - strength) > 1e-05 ) 
	{
		glLogger.debug("Different positions");
		return (false);	
	}
	else 
	{
		return (CxPosition::operator==(rhs));
	}
				
}

CxImpurity::~CxImpurity()
{
}

float CxImpurity::getPotential(const CxPosition &aufPunkt) const
{
	float retVal = 0.0f;
	float xP = aufPunkt.getXPosition();
	float yP = aufPunkt.getYPosition();
	float zaehler1 = (xP-getXPosition()) / getSigmaX();
	zaehler1 = zaehler1*zaehler1;
	float zaehler2 = (yP-getYPosition()) / getSigmaY();
	zaehler2 = zaehler2*zaehler2;
	float zaehler = zaehler1+ zaehler2;
	retVal = strength*exp(-1.0f*fabs(zaehler));
	return (retVal);
}

