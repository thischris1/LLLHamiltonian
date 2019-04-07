#include <LLLhamGaussianArray.h>
#include <utils/logger.hpp>
#include <ERRORS.h>

LLLhamGaussianArray::LLLhamGaussianArray(yosBasis* new_basis,int new_type,
										double bli_new, double a_new, 
										double b_new,
										double offset):
LLLhamiltonian(new_basis, new_type, bli_new, a_new, b_new),
m_impurities(impVector()),
m_offSet(offset),
xKmax(100),
yKmax(100),
m_tempArray(0)
{
	int Nm = new_basis->getNm();
	m_tempArray = new CMDCOMPLEX* [Nm];
	// Allocate the array 
	for (int index = 0; index < Nm; ++index) {
		CMDCOMPLEX *line = new CMDCOMPLEX[Nm];
		m_tempArray[index] = line;
	} 
	
	if (!preCalculateMatEls())
	{
		throw CxErrors(__FILE__,__LINE__);
	}
		
}


LLLhamGaussianArray::LLLhamGaussianArray(yosBasis* new_basis,int new_type,
										double bli_new, double a_new, 
										double b_new,
										double offset,
										impVector n_imps):
LLLhamiltonian(new_basis, new_type, bli_new, a_new, b_new),
m_impurities(n_imps),
m_offSet(offset),
xKmax(100),
yKmax(100),
m_tempArray(0)
{
	logImpurities();
	int Nm = new_basis->getNm();
	m_tempArray = new CMDCOMPLEX* [Nm];
	// Allocate the array 
	for (int index = 0; index < Nm; ++index) {
		CMDCOMPLEX *line = new CMDCOMPLEX[Nm];
		m_tempArray[index] = line;
	} 
	
	if (!preCalculateMatEls())
	{
		throw CxErrors(__FILE__,__LINE__);
	}
		
}

LLLhamGaussianArray::~LLLhamGaussianArray()
{
	if (m_tempArray)
	{
		for (int index = 0; index < Nm; ++index) {
			delete [] m_tempArray[index];
		} 
		delete [] m_tempArray;
	}
}

int LLLhamGaussianArray::matEl_perturb_nonzero(int i_index, int j_index, std::complex<double> *  retVal)
{
	int nonzero = 0;
  
  	const yosState *state1= (*LLL_basis)[i_index];
  	const yosState *state2=(*LLL_basis)[j_index];
	glLogger.debug("hamGaussianArray perturb nonzero state (%d), state (%d)", i_index, j_index);
   
  	CMDCOMPLEX ME;
  
  	if (i_index == j_index) // diagonal
    {
    
      nonzero = 1;
      ME += m_offSet;

      
      for (int i = 0; i < Ne; i++) 
		{
	  	/*!
	    	\todo here call the method which now calculates those values beforehand
	  	*/
	  

			ME = ME + getSPMatrixElement(state1->j(i), state1->j(i),0.0, 0.0);
			glLogger.debug("HamGaussianArray diagonal Matrix element is (%f + i(%f)", ME.real(), ME.imag());
		}
      
      
    }
  
  else  // off-diagonal 
    { 
      // find positions, where states differ
      int ind1=0, ind2=0;
      int diff1 = 0;           // counter for different WFs in state
      int diff2 = 0;
      int *diffs1 = new int [Ne];   // array that contains indices of WFs in state1 that are not contained in state2
      int *diffs2 = new int [Ne];   // ... vice versa ...
      
      while (ind1 < Ne && ind2 < Ne)
	{
	  if (state1->s1_nr(ind1) == state2->s2_nr(ind2)) {
	    ind1++; 
	    ind2++;	    
	  }
	  else {
	    if (state1->s1_nr(ind1) < state2->s2_nr(ind2)) 
	      {
		// all WF from state1 that are below st2.j(ind2) are not contained in state2
		while (state1->s1_nr(ind1) < state2->s2_nr(ind2) && ind1 < Ne) diffs1[diff1++] = ind1++;
	      }	      
	    else
	      {
		// all WF that are below st1.j(ind1) are not contained in state2
		while (state1->s1_nr(ind1) > state2->s2_nr(ind2) && ind2 < Ne) diffs2[diff2++] = ind2++;
	      }	    
	  } //of else if (st1.j(ind1) == st2.j(ind2))
	  
	  // if end of state1 reached -> rest of state2 is not contained in state1
	  if (ind1 == Ne)
	    while (ind2 < Ne) diffs2[diff2++] = ind2++;

	  // if end of state2 reached -> rest of state1 is not contained in state2
	  if (ind2 == Ne) 
	    while (ind1 < Ne) diffs1[diff1++] = ind1++;
	}
      
      // here diffs1/2 contains indices, that refer to positions in
      // state 1/2 have no corresponding wf in state 2/1

      
      if (diff1 == 1)
		{
		  nonzero = 1;
	  // parity = +-1
		  double parity = 1 - 2*((diffs1[0]+diffs2[0]) % 2);

	  
	  // hole
	  	int i = state1->j( diffs1[0] );
	  	int j = state2->j( diffs2[0] );
	  

	  /*!
	    \todo Call the function calculating m_pdHole... here
	  */

	  /* 
	     cm 26112004 dont need += ?
	  */
	  /*
	  ME = parity * (  m_pdHoleXeven[i+j] * (*m_pHoleYeven)( i, j )
			      + m_pdHoleXodd[i+j] * (*m_pHoleYodd)( i, j ) );
	  */
	  //ME = parity * totalMatEl(i, j, 0.0, 0.0); //! \todo alpha1 and alpha2 here!
		ME = parity*getSPMatrixElement(i,j, 0.0, 0.0);
		glLogger.debug("ME in hamGaussianArray offdiagonal PArt (%f + i %f)", ME.real(), ME.imag());

		}

      delete[] diffs1;
      delete[] diffs2;

    } // else  off diagonal
  
  	if (i_index == j_index)
    {
      if (fabs (ME.imag()) > 1e-20)
		{
		  	glLogger.error ("matrix-element on diagonal has imaginary part > 1e-20");
	  		glLogger.info ("ME(%d,%d) = (%3.30f,%3.30f)", i_index, j_index, ME.real(), ME.imag());
		}
      	ME = std::complex<double> (ME.real(), 0.0);
    }

  glLogger.debug("LLLhamGaussianArray::matEl_perturb_nonzero returns (%f), (%f)",ME.real(), ME.imag()),
  *retVal = ME;  
  
  
  return nonzero; // This matrix element is non-zero
	
}

bool LLLhamGaussianArray::clearImpurities()
{
	m_impurities.clear();
	return (preCalculateMatEls());
}
bool LLLhamGaussianArray::addImpurity(const CxImpurity &n_imp)
{
	
	m_impurities.push_back(n_imp);
	return (preCalculateMatEls());
	
}


/*!  \fn 
 * 
 * 
 */
CMDCOMPLEX LLLhamGaussianArray::getDiagonalElement(const int & index) const
{
	CMDCOMPLEX retVal;
	
return (retVal);	
}

CMDCOMPLEX LLLhamGaussianArray::getSPMatrixElement(const int &iIndex, const int &jIndex, double alpha1,double  alpha2) const
{
	glLogger.debug("Entering LLLhamGaussianArray::getSPMatrixElement (%d), (%d)", iIndex, jIndex);
	CMDCOMPLEX retVal = m_tempArray[iIndex][jIndex];
	
  /*
std::complex <double> 
  yEvenTmp,
  yOddTmp;

double 
  xEvenTmp =0.0,
  xOddTmp = 0.0;
   
 iPlusj(iIndex+jIndex,
	alpha2,
	xEvenTmp, 
	xOddTmp);
glLogger.debug("hamGaussianArray: SPMatrixElement: ipluusj returns xEven (%f)",xEvenTmp);
glLogger.debug("hamGaussianArray: SPMatrixElement: iplusj returns xOdd (%f)",xOddTmp);
 int iDiff= abs(iIndex-jIndex);
 iMinusj( yEvenTmp, yOddTmp, alpha1, iDiff);
 if (iIndex <  jIndex)
   {
     yEvenTmp = conj(yEvenTmp);
     yOddTmp = conj(yOddTmp);
   }
 
 CMDCOMPLEX retVal = (xEvenTmp*yEvenTmp) + (xOddTmp * yOddTmp); 
*/
glLogger.debug("SPMatrixelement returns (%f +i %f)", retVal.real(), retVal.imag());
 return (retVal);	
	
}


CMDCOMPLEX LLLhamGaussianArray::getSPAMAtrixElement(const int &iIndex, const int &jIndex, double alpha1,double  alpha2) const
{
	glLogger.debug("Entering LLLhamGaussianArray::getSPAMatrixElement (%d), (%d)", iIndex, jIndex);
	CMDCOMPLEX retVal; 
	
  
std::complex <double> 
  yEvenTmp,
  yOddTmp;

double 
  xEvenTmp =0.0,
  xOddTmp = 0.0;
   
 iPlusj(iIndex+jIndex,
	alpha2,
	xEvenTmp, 
	xOddTmp);
glLogger.debug("hamGaussianArray: SPAMatrixElement: ipluusj returns xEven (%f)",xEvenTmp);
glLogger.debug("hamGaussianArray: SPAMatrixElement: iplusj returns xOdd (%f)",xOddTmp);
 int iDiff= abs(iIndex-jIndex);
 iMinusj( yEvenTmp, yOddTmp, alpha1, iDiff);
 if (iIndex <  jIndex)
   {
     yEvenTmp = conj(yEvenTmp);
     yOddTmp = conj(yOddTmp);
   }
 
 retVal = (xEvenTmp*yEvenTmp) + (xOddTmp * yOddTmp); 

glLogger.debug("SPAMatrixelement returns (%f +i %f)", retVal.real(), retVal.imag());
 return (retVal);	
	
}

CMDCOMPLEX LLLhamGaussianArray::iMinusj(CMDCOMPLEX & MEYEvenPass, CMDCOMPLEX & MEYOddPass, double alpha1, int i_minus_j) const
{
	
	glLogger.debug("Entering GaussianArray iMinusj with (%d)",i_minus_j);
	CMDCOMPLEX retVal;
	std::complex<double> MEYeven(0.0, 0.0);
  	std::complex<double> MEYodd(0.0, 0.0);
	
	// loop over impurities
	for (unsigned int index = 0; index < m_impurities.size(); index++)
	{
		double strength = m_impurities[index].getStrength();
  		double sigmay = m_impurities[index].getSigmaY();
  		double xPos = m_impurities[index].getXPosition();
  		double yPos = m_impurities[index].getYPosition();

    	glLogger.info("Calculating impurity number (%d) with properties ", index);
   
  		glLogger.info( "x=(%f), y=(%f), sigmay=(%f)  Strength (%f)",
				 xPos,
				 yPos, 
				 sigmay,
				 strength);
		double 
    		alpha,
	    	beta,
    		a_to_b = a/b,
    		PhiEven = 0.0, 
    		PhiOdd  = 0.0;
  		alpha = 2.0 * M_PI * yPos;
  		beta =  M_PI * ( a_to_b / (2.0*Nm) + 
			   M_PI * sigmay* 
		   	sigmay);
		for (int l = -yKmax; l <= yKmax; l++)
    	{
      		PhiEven = fmod ( alpha * (-i_minus_j + 2*l*Nm) - alpha1 * 2*l, 2.0*M_PI );	  
		    MEYeven = MEYeven + 
				std::complex<double>(cos(PhiEven), sin(PhiEven)) * 
				exp( -beta * (i_minus_j - 2*l*Nm)*(i_minus_j - 2*l*Nm) )*
				fabs(strength);
    	} // end of kx loop 
//    	glLogger.debug("MEYEven after kx loop (%f), PhiEven= (%f +i%f)",PhiEven, MEYeven.real(), MEYeven.imag());

		for (int l = -yKmax; l < yKmax; l++)
    	{	  
      
		      PhiOdd = fmod ( alpha * (-i_minus_j + (2*l+1)*Nm) - alpha1 * (2*l+1), 2.0*M_PI );
		      MEYodd = MEYodd +
					(std::complex<double>(cos(PhiOdd), sin(PhiOdd)) * 
					 exp( -beta * (i_minus_j - (2*l+1)*Nm)*(i_minus_j - (2*l+1)*Nm) )*
					 fabs(strength));
    	} // end of ky loop
//    	glLogger.debug("after both loops (x and y) odd (%f +i  %f), even (%f +i%f) ", MEYodd.real(), MEYodd.imag(), MEYeven.real(), MEYeven.imag()); 
	  
	}
//	glLogger.debug("after loop over imps (%f), (%f) even (%f) (%f)", MEYodd.real(), MEYodd.imag(), MEYeven.real(), MEYeven.imag()); 
  	MEYEvenPass = MEYeven;
  	MEYOddPass = MEYodd;
  	glLogger.debug("Imunsj returns (%f), (%f)", MEYEvenPass.real(), MEYOddPass.real());
  	
	return (retVal);
	
}

MDOUBLE LLLhamGaussianArray::iPlusj(int  i_plus_j, double  alpha2, double &xEvenPass, double & yEvenPass) const
{
	glLogger.debug("LLLhamGaussianArray:: Entering iPlusj (%d)", i_plus_j);
	MDOUBLE retVal = 0.0;
	double 
		strength =0.0,
		sigmax= 0.0,
		sigmay = 0.0,
		xPos = 0.0;
	double a_to_b = a/b;
  	double eps = 2.0 * M_PI * Nm * a_to_b;
  	MDOUBLE  MEXeventmp = 0.0;
  	MDOUBLE  MEXoddtmp = 0.0;
  	double 
    	MEXeven =0.0,
    	MEXodd = 0.0;
	for ( unsigned int impIdx = 0; impIdx < m_impurities.size(); impIdx++)
	{ 
		
		strength = m_impurities[impIdx].getStrength();
  		sigmax = m_impurities[impIdx].getSigmaX();
	  	sigmay = m_impurities[impIdx].getSigmaY();
  		xPos = m_impurities[impIdx].getXPosition();
  		double sigma;
  		double prefact;
	  	glLogger.debug("Processing imp. no (%d) with x =(%f), strength=(%f)", impIdx, xPos,strength);
  		double imp2 = sigmax*sigmax;
	 
  		sigma = 1.0/(imp2 + 1.0/eps);
  		prefact = strength*sigmax*sigmay * sqrt(M_PI*sigma);
  		for (int r = -xKmax; r <= xKmax; r++)
    	{
      	/*
			Loop over different k's in the Yoshioka Wfct. sum
      	*/
	      	double argExp = r -xPos
				+ (i_plus_j+alpha2/M_PI)/(2.*Nm);
	      	MEXeventmp += exp ( -sigma * argExp*argExp );
	      	MEXoddtmp += exp ( -sigma * (argExp+0.5)*(argExp+0.5) );
     	}
   		MEXeven = MEXeven+MEXeventmp*prefact;
   		MEXodd = MEXodd+prefact*MEXoddtmp;
	}
   xEvenPass = MEXeven;
   yEvenPass = MEXodd;
		
	glLogger.debug("hamGaussianArray:iplusj, returns (%f),(%f)", xEvenPass, yEvenPass);	
	
	
	return (retVal);
}
bool LLLhamGaussianArray::logImpurities(void) const
{
	for (unsigned int impIdx = 0; impIdx < m_impurities.size(); ++impIdx) 
	{
		
	glLogger.info("	Processing imp. no (%d) with x =(%f), y =(%f), wx=(%f), wy =(%f), strength=(%f)", 
		impIdx, 
		m_impurities[impIdx].getXPosition(), 
		m_impurities[impIdx].getYPosition(),
		m_impurities[impIdx].getSigmaX(),
		m_impurities[impIdx].getSigmaY(),
		m_impurities[impIdx].getStrength() );
	}
return (true);	
}

bool LLLhamGaussianArray::preCalculateMatEls()
{
	glLogger.error("Start preCalculateMatEls()");
	for (int xVar = 0; xVar < Nm; ++xVar) {
		for (int yVar = 0; yVar < Nm; ++yVar) {
			m_tempArray[xVar][yVar] = getSPAMAtrixElement(xVar,yVar,0.0,0.0);
		}
		
	}
	return (true);
}
