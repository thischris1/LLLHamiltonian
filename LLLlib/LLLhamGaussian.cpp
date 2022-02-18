/*!
  \file LLLhamGaussian.cpp
  \brief holds the implementation of class LLLhamGaussian
  \author jakob sachs und christian mueller
*/


#include <LLLhamGaussian.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>

#include <Types.hpp>
/*!
  \brief constrcuts Hamiltonian for a single impurity system. Should go when teh other methods arte fully implementd
 \param new_basis the basis for the system
 \param new_type 
 \param bli_new correction factor due to finite width of the 2DEG
  \param a_new system size 
  \param b_new system size
  \param  xWidth Impurity parameters
  \param  xPosition
  \param  Strength
  \param yWidth
  \param yPosition
  \param Offset Energy offset
  \param xKmaxNew sum over k for the sp-wavefunctions x-direction
  \param yKmaxNew s.a. y-direction (ky)
 

 
 */

  
LLLhamGaussian::LLLhamGaussian( yosBasis *new_basis,
				int new_type,
				double bli_new,
				double a_new, 
				double b_new,
				double xWidth, 
				double xPosition,
				double Strength,
				double yWidth,
				double yPosition,
				double Offset,
				int xKmaxNew,
				int yKmaxNew ) 
  : LLLhamiltonian(new_basis,new_type,bli_new,a_new,b_new),
    m_dOffset(Offset),
	m_pdHoleXeven(0),
    m_pdHoleXodd(0),
    m_pHoleYeven(0) ,
    m_pHoleYodd(0),
    xKmax(xKmaxNew),
    yKmax(yKmaxNew),
    m_impurity(CxImpurity(xPosition, yPosition, xWidth, yWidth, Strength))
{
    glLogger.debug("Entering Constructor of hamGaussian with 1 impurity");
    calc_Hole ( 0, 0);
}

LLLhamGaussian::LLLhamGaussian(yosBasis *basis,  int type, double bli,
	       double a, double b, const CxImpurity &  newImpurities,
	        double offset, int xkmax, int ykmax ):
	        LLLhamiltonian(basis, type,bli,a,b),
    m_dOffset(offset),
	        m_pdHoleXeven(0),
m_pdHoleXodd(0),
    m_pHoleYeven(0) ,
    m_pHoleYodd(0),
    xKmax(xkmax),
    yKmax(ykmax),
    m_impurity(newImpurities)
{
	glLogger.debug("Entering Constructor of hamGaussian with type impurity");
    calc_Hole ( 0, 0);
}
/*!
  \brief Constructs Hamiltonian for a system with numOfImpurities Impurities.
  Parameter description as above
  \param  *basis the basis 
   \param  type Type of basis
   \param bli
   \param a size of system x Size
   \param  b size of system y size
   \param  newImpurities array of impurities
   \param numOfImpurities size of newImpurities / total number of impurities
   \param offset
   \param  xkmax
   \param  ykmax  

 */
 
/*!
  \brief Destructor

 */
LLLhamGaussian::~LLLhamGaussian() {

  delete m_pHoleYeven;
  delete m_pHoleYodd;
  delete [] m_pdHoleXeven;
  delete [] m_pdHoleXodd;
}
/*!
  \fn int LLLhamGaussian::matEl_perturb_nonzero(int i_st1,int i_st2, std::complex<double> *matEl)
  \param matEl returns the matrix element
  \param i_st1 Index of first State
  \param i_st2 Index of second state
*/

int LLLhamGaussian::matEl_perturb_nonzero(int i_st1,int i_st2, std::complex<double> *matEl)
{
  glLogger.debug("Entering  LLLhamGaussian::matEl_perturb_nonzero");
  int nonzero = 0;  
  CMDCOMPLEX ME(0.0,0.0);

  
  const yosState *state1= (*LLL_basis)[i_st1];
  const yosState *state2=(*LLL_basis)[i_st2];

   
  
  
  if (i_st1 == i_st2) // diagonal
    {
    
      nonzero = 1;
      ME += m_dOffset;

      for (int i = 0; i < Ne; i++) 
	{
	  /*!
	    \todo here call the method which now calculates those values beforehand
	  */
	  //  ME += m_pdHoleXEl[st1->j(i) + st1->j(i)] * (*m_pHoleYEl)(st1->j(i), st1->j(i));
	  
	  /*
	   * 	  	  ME +=  m_pdHoleXeven[2*state1->j(i)] * (*m_pHoleYeven)( state1->j(i), state1->j(i) )
	    + m_pdHoleXodd[2*state1->j(i)] * (*m_pHoleYodd)( state1->j(i), state1->j(i) );
	  	  ME = ME + totalMatEl(st1->j(i), st1->j(i), 0.0, 0.0); //! \todo alpha1 and alpha 2 here!
		*/
		ME = ME + totalMatEl(state1->j(i), state1->j(i), 0.0, 0.0);
		glLogger.debug(" HamGaussian matEl_perturb index (%d) = (%f + i%f)", i, ME.real(), ME.imag());
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
	    ind1++; ind2++;	    
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
	  	ME = parity * totalMatEl(i, j, 0.0, 0.0); //! \todo alpha1 and alpha2 here!
		}

    delete[] diffs1;
    delete[] diffs2;

    } // else  off diagonal
  
  if (i_st1 == i_st2)
    {
      if (fabs (ME.imag()) > 1e-20)
		{
	  	glLogger.error ("matrix-element on diagonal has imaginary part > 1e-20");
	  	glLogger.info ("ME(%d,%d) = (%3.30f,%3.30f)", i_st1, i_st2, ME.real(), ME.imag());
		}
      	ME = std::complex<double> (ME.real(), 0.0);
    }

  glLogger.debug("LLLhamGaussian::matEl_perturb_nonzero returns (%f), (%f)",ME.real(), ME.imag()),
  *matEl = ME;  
  // test *matEl = std::complex<double> (ME.real(), 0.0);
  
  return nonzero; // This matrix element is non-zero
}

  /////////////////////
 // private helpers //
/////////////////////


/*! 
  \fn void LLLhamGaussian::calc_Hole (double alpha1, double alpha2)
  \param alpha1 flux of solenoid 1 in units of h/e
  \param alpha2 flux of solenoid 2 in units of h/e
\brief precalculates factors for the impurity matrix
*/

void LLLhamGaussian::calc_Hole (double alpha1,double alpha2)
		
{
  glLogger.info ("LLLhamGaussian: calculating matrix-elements for hole\n");

  /*
    Set up local variables
   
  */

  //  double a_to_b = a/b;
  //  double eps = 2.0 * M_PI * Nm * a_to_b;
  double MEXeven = 0.0;
  double MEXodd = 0.0;
  if (m_pdHoleXeven == 0)
    {
      m_pdHoleXeven = new double[2*Nm-1];
      m_pdHoleXodd = new double[2*Nm-1];
    }

	  glLogger.info("HamGaussian Calculating impurity with properties X=(%f), y=(%f), sigmax (%f), (%f), Strength (%f)",
			 m_impurity.getXPosition(),
			 m_impurity.getYPosition(),
			m_impurity.getSigmaX(),
			m_impurity.getSigmaY(),
			 m_impurity.getStrength());
	 

  // calculation of part depending on i+j
  for (int i_p_j = 0; i_p_j < 2*Nm-1; i_p_j++)
    {
      MEXeven = 0.0;
      MEXodd = 0.0;
        double 
	    xPass = 0.0,
	    yPass = 0.0;
	
	  iplusj(i_p_j, alpha2,xPass, yPass);

      m_pdHoleXeven [i_p_j] = xPass; // MEXeven;
      m_pdHoleXodd [i_p_j] =  yPass; // MEXodd;
 
}
  
  //  double PhiEven, PhiOdd;

  std::complex<double> MEYeven;
  std::complex<double> MEYodd;
  if ( m_pHoleYeven == 0)
    {
      m_pHoleYeven = new BarrYEl(Nm);
      m_pHoleYodd = new BarrYEl(Nm);
    }

  // calculation of part depending on i-j
  for (int i_m_j = 0; i_m_j < Nm; i_m_j++)
    {
      	MEYeven = std::complex<double>(0.0, 0.0);
      	MEYodd = std::complex<double>(0.0, 0.0);
	  	iminusj(MEYeven, MEYodd, alpha1, i_m_j);
	 
	 
       	m_pHoleYeven->setAt (i_m_j, MEYeven);
       	m_pHoleYodd->setAt (i_m_j, MEYodd);
       
     	glLogger.debug(" MEYeven[ %d ] = (%3.12f, %3.15f)", i_m_j, real(MEYeven), imag(MEYeven));
	   	glLogger.debug (" MEYodd[ %d ] = (%3.12f, %3.15f)", i_m_j, real(MEYodd), imag(MEYodd));
	   
	 }

  glLogger.info ("...done calculating matrix elements of hole.");
}

// prints matrix elements on screen
void LLLhamGaussian::dumpElements ()
{
  for (int i=0; i<Nm; i++)
    for (int j=0; j<=i; j++)
      {

	std::complex<double> tmpV1 = 
	  m_pdHoleXeven[i+j]*(*m_pHoleYeven)(i,j) + m_pdHoleXodd[i+j]*(*m_pHoleYodd)(i,j);

  	glLogger.info ("<%d|V_hole|%d> = (%f, %f)", i, j, real(tmpV1), imag(tmpV1));
				   	
      }
}




/*!
  \param difference the difference between
  \brief calculation of part depending on i+j
*/
bool LLLhamGaussian::sumMatEl(int i_p_j, double & Xeven, double & Xodd) const
{

  double 
    alpha2 = 0.0;
  double a_to_b =  a/b;
  double eps = 2.0 * M_PI * Nm * a_to_b;
  double MEXeven = 0.0;
  double MEXodd = 0.0;
  
  
  
  MEXeven = 0.0;
  MEXodd = 0.0;
  
  glLogger.info("Calculating impurity with properties (%f), (%f), Strength (%f)",
		     m_impurity.getXPosition(),
		     m_impurity.getYPosition(),
		     m_impurity.getStrength());
      double sigma;
      double prefact;
      double  MEXeventmp = 0.0;
      double  MEXoddtmp = 0.0;
      double imp2 = 
	m_impurity.getSigmaX()*
	m_impurity.getSigmaX();
      
      sigma = 1.0/(imp2 + 1.0/eps);
      prefact =  m_impurity.getStrength()*
	m_impurity.getSigmaX()*
	m_impurity.getSigmaY()*
	sqrt(M_PI*sigma);
      for (int r = -xKmax; r <= xKmax; r++)
	{
	  /*
	    Loop over different k's in the Yoshioka Wfct. sum
	  */
	  double argExp = r -m_impurity.getXPosition()
	    + (i_p_j+alpha2/M_PI)/(2.*Nm);
	  MEXeventmp += exp ( -sigma * argExp*argExp );
	  MEXoddtmp += exp ( -sigma * (argExp+0.5)*(argExp+0.5) );
	}
      Xeven = MEXeven+MEXeventmp*prefact;
      Xodd = MEXodd+prefact*MEXoddtmp;
  
      return (true);
      
    }

/*!


*/
bool LLLhamGaussian::diffMatEl(int i_m_j, complex <double>& YEven, complex <double>& YOdd) const
{

  double PhiEven, PhiOdd;
  double a_to_b = a/b;
  double
    alpha1 = 0.0;
  //  double eps = 2.0 * M_PI * Nm * a_to_b;
  std::complex<double> MEYeven = std::complex<double>(0.0, 0.0);
  std::complex<double> MEYodd  = std::complex<double>(0.0, 0.0);

       /*
	Calculating a few prefactors from the impurities properties
      */
      glLogger.info("Calculating impurity with properties (%f), (%f), Strength (%f)",
		     m_impurity.getXPosition(),
		     m_impurity.getYPosition(),
		     m_impurity.getStrength());
      double 
	alpha,
	beta;
      alpha = 2.0 * M_PI * m_impurity.getYPosition();
      beta =  M_PI * ( a_to_b / (2.0*Nm) + 
		       M_PI * m_impurity.getSigmaY()* 
		       m_impurity.getSigmaY());
      for (int l = -yKmax; l <= yKmax; l++)
		{
	  		PhiEven = fmod ( alpha * (-i_m_j + 2*l*Nm)- alpha1 * 2*l, 2.0*M_PI );	  
	  		MEYeven = MEYeven + 
	    		CMDCOMPLEX(cos(PhiEven), sin(PhiEven)) * 
	    		exp( -beta * (i_m_j - 2*l*Nm)*(i_m_j - 2*l*Nm) )*
	    		fabs(m_impurity.getStrength());
		}
      
      for (int l = -yKmax; l < yKmax; l++)
		{	  
	  
	  		PhiOdd = fmod ( alpha * (-i_m_j + (2*l+1)*Nm), 2.0*M_PI  );
	  		MEYodd = MEYodd +
	    	(CMDCOMPLEX(cos(PhiOdd), sin(PhiOdd)) * 
	     	exp( -beta * (i_m_j - (2*l+1)*Nm)*(i_m_j - (2*l+1)*Nm) )*
	     	fabs(m_impurity.getStrength()));
		}
    
  
 	 glLogger.debug(" MEYeven[ %d ] = (%3.12f, %3.15f)", i_m_j, real(MEYeven), imag(MEYeven));
  	glLogger.debug (" MEYodd[ %d ] = (%3.12f, %3.15f)", i_m_j, real(MEYodd), imag(MEYodd));
  	YEven = MEYeven;
  	YOdd = MEYodd;
  	return (true);
  
}




double LLLhamGaussian::iplusj(int  i_plus_j,
			      double  alpha2,
			      double &xEvenPass,
			      double &yEvenPass) const
{
//	glLogger.debug("LLLhamGaussin:: Entering iplusj with (%d)", i_plus_j);
  double strength = m_impurity.getStrength();
  double sigmax = m_impurity.getSigmaX();
  double sigmay = m_impurity.getSigmaY();
  double xPos = m_impurity.getXPosition();

 double a_to_b = a/b;
  double eps = 2.0 * M_PI * Nm * a_to_b;
  double retVal = 0.0;
  double sigma;
  double prefact;
  double  MEXeventmp = 0.0;
  double  MEXoddtmp = 0.0;
  double 
    MEXeven =0.0,
    MEXodd = 0.0;
  double imp2 = 
    sigmax*sigmax;
	 
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
  
   xEvenPass = MEXeven;
   yEvenPass = MEXodd;
 
 	glLogger.debug("HamGaussion iplusj = (%d) (%f), (%f)", i_plus_j, xEvenPass, yEvenPass);
   return (retVal);
}


bool LLLhamGaussian::iminusj(  
			     std::complex <double> & MEYEvenPass,
			     std::complex <double> & MEYOddPass,
			     double alpha1,	       
			     int i_minus_j) const
{
//	glLogger.debug("LLLhamGaussian::iminusj Entering with (%d)", i_minus_j); 
	double strength = m_impurity.getStrength();
  	double sigmax = m_impurity.getSigmaX();
  	double sigmay = m_impurity.getSigmaY();
  	double xPos = m_impurity.getXPosition();
  	double yPos = m_impurity.getYPosition();

  std::complex<double> MEYeven(0.0, 0.0);
  std::complex<double> MEYodd(0.0, 0.0);
  glLogger.debug("LLLhamGauss: iminusj Calculating impurity with properties :x=(%f), y=(%f) sigmax =(%f), sigmay=(%f) , Strength (%f)",
		 xPos,
		 yPos,
		 sigmax,
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
    }

  for (int l = -yKmax; l < yKmax; l++)
    {	  
      
      PhiOdd = fmod ( alpha * (-i_minus_j + (2*l+1)*Nm) - alpha1 * (2*l+1), 2.0*M_PI );
      MEYodd = MEYodd +
	(std::complex<double>(cos(PhiOdd), sin(PhiOdd)) * 
	 exp( -beta * (i_minus_j - (2*l+1)*Nm)*(i_minus_j - (2*l+1)*Nm) )*
	 fabs(strength));
    }
  MEYEvenPass = MEYeven;
  MEYOddPass = MEYodd;

  
//glLogger.debug("Imunsj returns (%f), (%f)", MEYEvenPass.real(), MEYOddPass.real());
  return (true);
}


 std::complex <double> LLLhamGaussian::singleME ( int  iIndex,
				  int  jIndex,
				  double alpha1,
				  double alpha2) const
{
glLogger.debug("hamGaussian singleME enter with (%d), (%d)", iIndex, jIndex);   
std::complex <double> 
  yEvenTmp,
  yOddTmp,
  retVal;

double 
  xEvenTmp =0.0,
  xOddTmp = 0.0;
 int  plus = iIndex + jIndex; 
 iplusj(plus,
	alpha2,
	xEvenTmp, 
	xOddTmp);

 int iDiff = abs(iIndex-jIndex);
 iminusj( yEvenTmp, yOddTmp,
	 alpha1,
	 iDiff);
 if (iIndex <  jIndex)
   {
     yEvenTmp = conj(yEvenTmp);
     yOddTmp = conj(yOddTmp);
   }
 
 retVal = (xEvenTmp*yEvenTmp) + (xOddTmp * yOddTmp); 
glLogger.debug("Hamgaussian singleME returns (%f +i %f)", retVal.real(), retVal.imag());
 return (retVal);
}



std::complex <double> LLLhamGaussian::totalMatEl(int iIndex,
						 int jIndex,
						 double alpha1,
						 double alpha2) const
{
  
  CMDCOMPLEX retVal = singleME(iIndex, jIndex, alpha1, alpha2);
  glLogger.debug("Hamgaussian totalMatel returns (%f +i %f)", retVal.real(), retVal.imag());
  return (retVal);

}

bool LLLhamGaussian::setImpurity(const CxImpurity& rhs)
{
	m_impurity=rhs;
	return (true);
}
