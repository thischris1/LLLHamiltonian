/*!
  \file LLLhamBarrier.cpp
  \author Moritz helias
  \brief Holds the implementation of class LLLhamBarrier
  
*/


#include <LLLhamBarrier.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>

LLLhamBarrier::LLLhamBarrier( yosBasis *new_basis,
			      int new_type,
			      double bli_new,
			      double a_new, 
			      double b_new,
			      double Strength, 
			      double xWidth, 
			      double xPosition,
			      double holeStrength,
			      double holeWidth,
			      double holePosition,
			      double Offset,
			      double alpha1,
			      double alpha2,
			      int iPotType,
			      int xKmax,
			      int yKmax ) : LLLhamiltonian(new_basis,new_type,bli_new,a_new,b_new)
{
  m_iPotType = iPotType;
  m_dBarrStrength = Strength;
  m_dOffset = Offset;
  m_pdBarrXEl = NULL;
  m_pBarrYEl = NULL;
  m_pdHoleXeven = NULL;
  m_pdHoleXodd = NULL;
  m_pHoleYeven = NULL;
  m_pHoleYodd = NULL;

  m_iBarrX = (int) (xPosition*Nm+0.5);

  // need to pass 1.0 as height to use this matrixelement as a factor for
  // the delta shaped hole
  // will be multiplied in matElnonzero with actual height of barrier
  calc_BarrierX_el (xPosition, 1.0, xWidth, alpha2, xKmax);
  calc_Hole (holeStrength, xPosition, holePosition, xWidth, holeWidth, a_new/b_new, alpha1, alpha2, xKmax, yKmax);
  calc_BarrierY_el (holePosition, holeStrength, holeWidth, alpha1,  yKmax);

  //dumpElements (); 
}

LLLhamBarrier::~LLLhamBarrier() {
  delete m_pBarrYEl;
  delete [] m_pdBarrXEl;

  delete m_pHoleYeven;
  delete m_pHoleYodd;
  delete [] m_pdHoleXeven;
  delete [] m_pdHoleXodd;
}

/*!
  \fn int LLLhamBarrier::matEl_perturb_nonzero(int i_st1,int i_st2, std::complex<double> *matEl)
  \param matEl returns the matrix element
  \param i_st1 State to be used for calculation of matrix element
  \param i_st2
*/

int LLLhamBarrier::matEl_perturb_nonzero(int i_st1,int i_st2, std::complex<double> *matEl)
{
  glLogger.info("Entering  LLLhamBarrier::matEl_perturb_nonzero");
  int nonzero = 0;  
  std::complex<double> ME;


  st1=(*LLL_basis)[i_st1];
  st2=(*LLL_basis)[i_st2];

   
  ME = std::complex<double>(0.0,0.0);
  
  if (i_st1 == i_st2) // diagonal
    {

      // barrier parallel to y-axes only diagonal-terms
      switch (m_iPotType)
	{
	case 0: // gaussian
	  for(int i=0; i<Ne; i++) 
	    ME += m_dBarrStrength*m_pdBarrXEl [ st1->j(i) ];
	  break;
	  
	case 1: // delta
	  for(int i=0; i<Ne; i++) 
	    if (st1->j(i) == m_iBarrX) 
	      ME += m_dBarrStrength;
	  break;
	}

    nonzero = 1;

    ME += m_dOffset;

    // Barrier parallel to x-axes
    // for (int i = 0; i < Ne; i++) 
    //  MEY += (*m_pBarrYEl) ( st1->j(i), st1->j(i) );

    // hole
    switch (m_iPotType)
      {
      case 0:
	for (int i = 0; i < Ne; i++) 
	  //  ME += m_pdHoleXEl[st1->j(i) + st1->j(i)] * (*m_pHoleYEl)(st1->j(i), st1->j(i));
	  ME +=  m_pdHoleXeven[2*st1->j(i)] * (*m_pHoleYeven)( st1->j(i), st1->j(i) )
	       + m_pdHoleXodd[2*st1->j(i)] * (*m_pHoleYodd)( st1->j(i), st1->j(i) );
	break;

      case 1:
	for (int i = 0; i < Ne; i++) 
	  if (st1->j(i) == m_iBarrX) 
	    ME += (*m_pBarrYEl) ( st1->j(i), st1->j(i) );
	break;	
      }

    }

  else  // off-diagonal 
    { 
      // find positions, where states differ
      int ind1=0, ind2=0;
      int diff1 = 0;           // counter for different WFs in state
      int diff2 = 0;
//      int *diffs1 = new int [Ne];   // array that contains indices of WFs in state1 that are not contained in state2
//      int *diffs2 = new int [Ne];   // ... vice versa ...
	std::vector<int>  diffs1(Ne);
	std::vector<int>  diffs2(Ne);
      while (ind1 < Ne && ind2 < Ne)
	{
	  if (st1->s1_nr(ind1) == st2->s2_nr(ind2)) {
	    ind1++; ind2++;	    
	  }
	  else {
	    if (st1->s1_nr(ind1) < st2->s2_nr(ind2)) 
	      {
		// all WF from state1 that are below st2.j(ind2) are not contained in state2
		while (st1->s1_nr(ind1) < st2->s2_nr(ind2) && ind1 < Ne) diffs1[diff1++] = ind1++;
	      }	      
	    else
	      {
		// all WF that are below st1.j(ind1) are not contained in state2
		while (st1->s1_nr(ind1) > st2->s2_nr(ind2) && ind2 < Ne) diffs2[diff2++] = ind2++;
	      }	    
	  } //of else if (st1.j(ind1) == st2.j(ind2))

	  // if end of state1 reached -> rest of state2 is not contained in state1
	  if (ind1 == Ne)
//	    while (ind2 < Ne) diffs2[diff2++] = ind2++;
		while (ind2 < Ne) 
		{
			diffs2.at(diff2++) = ind2++;
		}
	  // if end of state2 reached -> rest of state1 is not contained in state2
	  if (ind2 == Ne) 
	    //while (ind1 < Ne) diffs1[diff1++] = ind1++;
		while (ind1 < Ne) 
		{
			diffs1.at(diff1++) = ind1++;
		}
	}
      
      // here diffs1/2 contains indices, that refer to positions in
      // state 1/2 have no corresponding wf in state 2/1

      
      if (diff1 == 1)
	{
	  nonzero = 1;
	  // parity = +-1
	  double parity = 1 - 2*((diffs1[0]+diffs2[0]) % 2);

	  // barrier parallel to x-axes
	  //MEY = parity * (*m_pBarrYEl) ( st1->j( diffs1[0] ), st2->j( diffs2[0] ) );

	  // hole
	  int i = st1->j( diffs1[0] );
	  int j = st2->j( diffs2[0] );

	  /*
	  std::cout << "state1: |";
	  for (int n=0; n<Ne; n++) std::cout << st1->j(n) << ' ';
	  std::cout << ">\n";

	  std::cout << "state2: |";
	  for (int n=0; n<Ne; n++) std::cout << st2->j(n) << ' ';
	  std::cout << ">\n";

	  std::cout << "differ at " << diffs1[0] << ' ' << diffs2[0] << " parity=" << parity << '\n';
	  */
	  switch (m_iPotType)
	    {
	    case 0:
	      ME += parity * (  m_pdHoleXeven[i+j] * (*m_pHoleYeven)( i, j )
			      + m_pdHoleXodd[i+j] * (*m_pHoleYodd)( i, j ) );
	      break;

	    case 1:
	      if (i == m_iBarrX || j == m_iBarrX)
		//ME += parity * 0.5 * ((*m_pBarrYEl) (i, j) + (*m_pBarrYEl) (j, i));
		ME += parity * (*m_pBarrYEl) (i, j);
		
	      break;
	    }


	}

     // delete[] diffs1;
      //delete[] diffs2;

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

  glLogger.info("LLLhamBarrier::matEl_perturb_nonzero returns (%f), (%f)",ME.real(), ME.imag()),
  *matEl = ME;  
  // test *matEl = std::complex<double> (ME.real(), 0.0);
  
  return nonzero; // This matrix element is non-zero
}

  /////////////////////
 // private helpers //
/////////////////////

  // evalutates one-particle matrix-elements of barrier along x-axes (parallel to y-axes)
  // p: Position in units of a
  // s: peak strength in units of e^2/(epsilon l0)
  // w: width in units of a
  // alpha2: Flux of solenoid2 in units of e/h

void LLLhamBarrier::calc_BarrierX_el (double p, double s, double w, double alpha2, int kmax) {
  double w2 = w*w;
  double tmp = 2.0*M_PI*Nm*(a/b);

  double Lambda = tmp/(1.0+tmp*w2);
  double Sigma = s/sqrt(1.0+1.0/(tmp*w2));

  double MEreal, exponent;
  
  // barrier parallel to y-axes only diagonal-terms
  glLogger.info ("calculating matrix-elements for BarrierX\n");

  m_pdBarrXEl = new double[Nm];

  
  for (int i=0; i<Nm; i++)
    {
      MEreal=0;
      for (int k=-kmax; k<=kmax; k++) 
	{
	  exponent = ( ( (double)i + alpha2/(2.0*M_PI) ) / Nm + k - p);
	  MEreal += exp( - Lambda * exponent*exponent );
	}
      m_pdBarrXEl[i] = Sigma*MEreal;
      //std::cout << i << ' ' << Sigma*MEreal << '\n';
    }
}


  // evalutates one-particle matrix-elements of barrier along y-axes
  // p: Position in units of b
  // s: peak strength in units of e^2/(epsilon l0)
  // w: width in units of b
  // alpha1: flux of solenoid #1 in units of h/e
  // kmax: maximum k to sum over
void LLLhamBarrier::calc_BarrierY_el (double p, double s, double w, double alpha1, int kmax) 
{
  double sigma = s*w*sqrt(M_PI);
  double lambda = M_PI * ((a/b)/(2.0*Nm) + M_PI*w*w);
  double peakheight;
  
  std::complex<double> ME;

  glLogger.info ("calculating matrix-elements for BarrierY\n");

  // First determine effective height of barrier
  // because due to the overlapp from neighboring cells the peak height may have to be
  // reduced to achieve the desired effective peak height.
  
  // calculate peak height of superposition
  // from all neighboring cells
  peakheight = 0.0;
  for (int k = -kmax; k < kmax; k++)
    {
      peakheight += exp ( - k*k/(w*w) );
    }

  glLogger.info("correction factor for peak height is %f \n", 1.0/peakheight);

  sigma = sigma / peakheight;

  // Nm elements to calculate because matrix elements only depend on i-j

  m_pBarrYEl = new BarrYEl(Nm);  


  for (int i_j = 0; i_j < Nm; i_j++)
    {
      
      // <i|V(y)|j>
      ME = std::complex<double>(0.0,0.0);

      for (int k=-kmax; k<=kmax; k++) {
	
	double arg = fmod ( 2.0*M_PI*( -i_j + k*Nm ) * p + alpha1*k, 2.0*M_PI ); // only fractional part
	//assert (arg >= -2*M_PI && arg <= 2*M_PI);

	ME += exp( -lambda * (k*Nm-i_j)*(k*Nm-i_j) ) * std::complex<double>( cos(arg), sin(arg) );      

	// for special case p = 0.5
	//double parity = 1 - ((i_j) % 2) * 2;
	//ME += exp( -lambda * (k*Nm-i_j)*(k*Nm-i_j) ) * parity;
      }

      if (ME.imag() != 0.0) 
        {
	  std::cout << ME << '\n';
        }

      m_pBarrYEl->setAt( i_j, sigma*ME );
    }

  glLogger.info ("...done calculating matrix elements BarrierY\n");

}
/*! 
  \fn void LLLhamBarrier::calc_Hole (double s, double px, double py, double wx,double wy, double a_to_b, double alpha1, double alpha2, int kxmax, int kymax)
  \param s Strength
  \param px position in units of a
  \param py position in units of b
  \param wx width in units of a
  \param wy width in units of b
  \param a_to_b = a/b aspect ratio
  \param alpha1 flux of solenoid 1 in units of h/e
  \param alpha2 flux of solenoid 2 in units of h/e
  \param kxmax maximuk k to sum over in x-sum
  \param kymax maximuk k to sum over in y-sum
*/

void LLLhamBarrier::calc_Hole (double s, 
			       double px, 
			       double py, 
			       double wx, 
			       double wy, 
			       double a_to_b, 
			       double alpha1,
			       double alpha2,
			       int kxmax, 
			       int kymax)
{
  glLogger.info ("calculating matrix-elements for hole\n");

    double eps = 2.0 * M_PI * Nm * a_to_b;
  double sigma = 1.0 / (wx*wx + 1.0/eps);
  double prefact = s * wx * wy * sqrt(M_PI*sigma);
  double MEXeven = 0.0;
  double MEXodd = 0.0;

  m_pdHoleXeven = new double[2*Nm-1];

  m_pdHoleXodd = new double[2*Nm-1];


  // calculation of part depending on i+j
  for (int i_p_j = 0; i_p_j < 2*Nm-1; i_p_j++)
    {
      MEXeven = 0.0;
      MEXodd = 0.0;
      for (int r = -kxmax; r <= kxmax; r++)
	{
	  /*
	    Loop over different k's in the Yoshioka Wfct. sum
	  */
	  double argExp = r - px + (i_p_j+alpha2/M_PI)/(2.*Nm);
	  MEXeven += exp ( -sigma * argExp*argExp );
	  MEXodd += exp ( -sigma * (argExp+0.5)*(argExp+0.5) );
	}
      
      m_pdHoleXeven [i_p_j] = prefact * MEXeven;
      m_pdHoleXodd [i_p_j] = prefact * MEXodd;
      glLogger.debug("MEXeven (%d) = (%f)",  
		     i_p_j, m_pdHoleXeven [i_p_j]);
      glLogger.debug("MEXodd (%d) = (%f)",
		     i_p_j, m_pdHoleXodd [i_p_j]);
    }
  


  double alpha = 2.0 * M_PI * py;
  double beta = M_PI * ( a_to_b / (2.0*Nm) + M_PI * wy*wy );

  double PhiEven, PhiOdd;

  std::complex<double> MEYeven;
  std::complex<double> MEYodd;

  m_pHoleYeven = new BarrYEl(Nm);
  m_pHoleYodd = new BarrYEl(Nm);


  // calculation of part depending on i-j
  for (int i_m_j = 0; i_m_j < Nm; i_m_j++)
    {
      MEYeven = std::complex<double>(0.0, 0.0);
      MEYodd = std::complex<double>(0.0, 0.0);

      for (int l = -kymax; l <= kymax; l++)
	{
	  PhiEven = fmod ( alpha * (-i_m_j + 2*l*Nm) - alpha1 * 2*l, 2.0*M_PI );	  
	  MEYeven += 
	    std::complex<double>(cos(PhiEven), sin(PhiEven)) * 
	    exp( -beta * (i_m_j - 2*l*Nm)*(i_m_j - 2*l*Nm) );

	}

      for (int l = -kymax; l < kymax; l++)
	{	  

	  PhiOdd = fmod ( alpha * (-i_m_j + (2*l+1)*Nm) - alpha1 * (2*l+1), 2.0*M_PI );
	  MEYodd += 
	    std::complex<double>(cos(PhiOdd), sin(PhiOdd)) * 
	    exp( -beta * (i_m_j - (2*l+1)*Nm)*(i_m_j - (2*l+1)*Nm) );
	}


      m_pHoleYeven->setAt (i_m_j, MEYeven);
      m_pHoleYodd->setAt (i_m_j, MEYodd);
      
      glLogger.debug(" MEYeven[ %d ] = (%3.12f, %3.15f)", i_m_j, real(MEYeven), imag(MEYeven));
      glLogger.debug (" MEYodd[ %d ] = (%3.12f, %3.15f)", i_m_j, real(MEYodd), imag(MEYodd));

    }

  glLogger.info ("...done calculating matrix elements of hole.");
}

// prints matrix elements on screen
void LLLhamBarrier::dumpElements ()
{
  for (int i=0; i<Nm; i++)
    for (int j=0; j<=i; j++)
      {

	std::complex<double> tmpV1 = 
	  m_pdHoleXeven[i+j]*(*m_pHoleYeven)(i,j) + m_pdHoleXodd[i+j]*(*m_pHoleYodd)(i,j);

  	glLogger.info ("<%d|V_hole|%d> = (%f, %f)", i, j, real(tmpV1), imag(tmpV1));
		

	std::complex<double> tmpV2 = 0.5*( m_pdBarrXEl[i]*(*m_pBarrYEl)(i,j) + m_pdBarrXEl[j]*(*m_pBarrYEl)(i,j) );

	glLogger.info ("<%d|V_x*V_y|%d> = (%f, %f)", i, j, real(tmpV2), imag(tmpV2));
		   	
      }
}
