  /*!
  \file yosState.cpp
  \brief Implementation of the class yosState
  \author Karel Vyborny Christian Mueller
*/
#include <yosState.hpp>

#include <utils/logger.hpp>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>


/*!
  \fn int equalSpat(const yosState *st1,const yosState *st2)
  \param st1
  \param st2
  (a friend function of yosState)
   Compares two yosState's st1, st2.
   \return 1 if the spatial parts are the same (i.e. all j's are the same) and the states
   are of the same type (see equalType()). Returns 0 otherwise.
*/
//! \todo This is not a friend functioN!!! It should go into yosState !!

int equalSpat(const yosState *st1,const yosState *st2)
{
  int n=st1->getNe();
  if(!equalType(st1,st2)) return 0;
  for(int i=0;i<n;i++) if(st1->j(i)!=st2->j(i)) return 0;
  return 1;
}


/* (a friend function of yosState)
   Compares two yosState's st1, st2. Returns 1 if Ne's are
   the same and Nm's are the same. Returns 0 otherwise.
*/
//! \todo This is not a friend functioN!!! It should go into yosState !! 

int equalType(const yosState *st1,const yosState *st2)
{
  if(st1->getNe()!=st2->getNe()) return 0;
  if(st1->getNm()!=st2->getNm()) return 0;
  return 1;
}

/*!  
  \fn yosState::yosState(int new_Ne,int new_Nm,int *new_j)
  \param new_Ne
  \param new_Nm
  \param new_j
  \brief Constructor for a new Yoshioka type state
*/
yosState::yosState(int new_Ne,int new_Nm,int *new_j):
Ne_(new_Ne),
Nm_(new_Nm),
spinYes(0),
j_(new int[Ne_]),
spin_(0),
aToB(1.0),
alpha1_(0.0),
alpha2_(0.0)
{
  
  
  for(int i=0;i<Ne_;i++) j_[i]=new_j[i];  
  }

/*!
\fn yosState::yosState(int new_Ne,int new_Nm,int *new_j,int *new_spin)
 \param new_Ne
  \param new_Nm
  \param new_j
  \param new_spin
  \brief Constructor for a new Yoshioka type state
*/
yosState::yosState(int new_Ne,int new_Nm,int *new_j,int *new_spin):
Ne_(new_Ne),
Nm_(new_Nm),
spinYes(1),
j_(new int[Ne_]),
spin_(0),
aToB(1.0),
alpha1_(0.0),
alpha2_(0.0)
{
  
  spin_=new int[Ne_];
  spinYes=1;
  for(int i=0;i<Ne_;i++) 
  	{
  		j_[i]=new_j[i];
  		spin_[i]=new_spin[i];
  	}
 
  }
/*!
  \fn yosState::yosState( yosState& rhs)
  \brief Copy constructor: deep. (originally intended as shallow; detach() 
  would have to be implemented!)
  \param rhs the State to be copied

  \todo detach() & shallow copy-ctor

To generate a deep copy use the detach() method (Idea stolen from qt)
*/

yosState::yosState(const yosState& rhs):
spin_(0)
{
	if (&rhs != this)
	{ 
  	Ne_ = rhs.getNe();
  	Nm_ = rhs.getNm();
  	alpha1_ = rhs.getSolenoidFlux1();
  	alpha2_ = rhs.getSolenoidFlux2();
  	aToB = rhs.getAspect();
  	spinYes = rhs.hasSpin();
  	/* This is a shallow copy 
  	j_ = rhs.getjStates();
  	if (spinYes == 1)
    {
      spin_ = rhs.getSpinStates();
    }
  // Reference counter 
  reference++;
  */
  
  /* This is a deep copy */
  	j_=new int[Ne_];
  	for(int i=0;i<Ne_;i++) j_[i]=rhs.getjStates()[i];
  	if (spinYes == 1)
    	{
      	spin_ = new int[Ne_];
      	for(int i=0;i<Ne_;i++) 
      		{
      			spin_[i]=rhs.getSpinStates()[i];
      		}
    	}
  	}
}

/*!
  \fn yosState::~yosState()
  \brief Destructor
*/
yosState::~yosState()
{
  
    
    delete[] j_;
    j_=0;
    delete[] spin_;
    spin_ =0;
     
}


/*!
  \fn yosState::operator = (const yosState &rhs)
  \brief Copying via copy-constructor.
*/

yosState & yosState::operator = (const yosState &rhs)
{
	if (&rhs != this)
	{
		
		
	}
  yosState *yos = new yosState(rhs); // call the copy constructor
  return *yos; 
}


/*!
  \fn const void yosState::print() const
  \brief prints the state in occupation number style | 1,2,3 , ... >
*/
void yosState::print() const
{
  printf("|");
  for(int i=0;i<Ne_;i++)
    {
      if(i>0) printf(" ");
      printf("%d",j_[i]);
      if(spinYes)
      {
    	  if(spin_[i]) printf("+,");
      }
	else
	  printf("-,");
    }
  printf(">");//printf("\n");
}
void yosState::fprint_x(FILE *f) const
{
  char *str = new char[Nm_+3];
  sPrint_occ(str,Nm_+3);
  fprintf(f,"%s",str);  
  delete []str;
}


void yosState::fprint(FILE *f) const
{
  if (f != 0)
    {
      fprintf(f,"|");
    }
  for(int i=0;i<Ne_;i++)
    {
      if(i>0) 
	{
	  if (f != 0)
	    {
	      fprintf(f," ");
	      fprintf(f,"%d",j_[i]);
	    }
	 
	  if(spinYes) 
	    {
	      if(spin_[i]) 
		{
		  fprintf(f,"+,");
		}
	      else
		{
		  fprintf(f,"-,");
		}
	    }
	}
    }
  fprintf(f,">");
}

/*!
  \fn bool yosState::compare (const yosState & rhs) const
  \brief compares two states (occ number wise
  \return true if equal
  \todo enhance that to find all states which are equal even if interchanged
  (beware of sign!)
*/
bool yosState::compare (const yosState & rhs) const
{

  if (rhs.getNe() != Ne_  ||
      rhs.getNm() != Nm_  ||
      (rhs.hasSpin() !=  spinYes ) )
    {
      glLogger.info( "Compare returns false, for obvious reasons");; 
      return false;
    }
    if (spinYes && !spin_)
    {
    	  return (false);
    }
   
  for (int iIdx = 0; iIdx < Ne_; iIdx++)
    {
      //      cerr << "Compare (" << rhs.j(iIdx) <<") and ("<< j_[iIdx] <<")\n";
      if (rhs.j(iIdx) != j_[iIdx])
	{
	  return false;
	}
      if ( (spinYes == 'y')
	   && (spin_[iIdx] != spin(iIdx)))
	{
//	  std::cerr << "Compare two spinstates";
	  return false;
	}
    }
      return true;
      
}

/*!

\fn 


*/

bool yosState::operator == (const yosState& rhs) const
{
  return (compare (rhs));
}



/*!
  \fn  bool yosState::getSingleDiffOccNumber ( const yosState &rhs, int& rIndex, int& lIndex,int& sign  ) const
\brief This methods compares two yosStates. It returns true if they differ in only one occ-number (and the rest of tha parameters are equal of course). Should be used in population of matrices (for example single particle operators).
\return true if the states differ in only one occupationsnumber, false otherwise see also the parameters for returned values
\param rhs the state to be compared with
\param rIndex the Index in the rhs State which is differing (on Exit)
\param  lIndex the Index in this state  which is differing (on Exit)
\param sign the sign (-1 or 1) depending on the numbers of neede cyclical permutations to put the states into proper order
\todo extend this to spinstates 
*/

bool yosState::getSingleDiffOccNumber ( const yosState &rhs,  
					int& rIndex,  
					int& lIndex,  
					int& sign  ) const
{

if (hasSpin())
    {
    glLogger.error( "getSingleDiffOccNumber cant handle spinstates yet. Leave");
      return false;
    }
	
 if (*this == rhs)
    { 
      return false;
    }
 
 if (rhs.getNe() != Ne_  ||
     rhs.getNm() != Nm_  ||
     (rhs.hasSpin() !=  spinYes ) )
   {
     glLogger.info( "Compare returns false, for obvious reasons\n"); 
     return false;
   }
 int 
   countState = 0,
   diffState = 0;
   

  std::vector <int>  takenRSide;
  std::vector <int> takenLSide;
 
 //int takenRSide[Ne_-1];
 // int takenLSide[Ne
 
 for (int iIdx =0; iIdx < Ne_; iIdx++)
   {
     bool found = false;
     for (int iIdy = 0; iIdy < Ne_; iIdy++)
       {
	 // find an equal element to j_[iIdx] 
	 if (j_[iIdx] == rhs.j(iIdy))
	   {
	     // Check if iIdy is already taken
	     bool taken = false;
	     for (int iTemp =0; iTemp < countState; iTemp++)
	       {
		 if (iIdy == takenRSide[iTemp])
		   {
		     taken = true;
		     break;
		   }
	       }
	     if (!taken)
	       {
		 // Copy the indices into lState and Rstate
		 takenLSide.push_back(iIdx);
		 takenRSide.push_back(iIdy);
		 found = true;
		 break; // Get to next element
	       }
	   }
       }
     if (!found) 
       {
	 diffState++;
	 if (diffState > 1)
	   {
	     return false;
	   }
       }
   }
 // If we get here,  states are exactly diferrent at 1 occ number. The indices for the pairs are in takenLSide and takenRSide. Now find the missing indices
 
     
 for (int iTemp =0; iTemp < Ne_; iTemp++)
   {
     std::vector<int>::iterator it = takenRSide.begin();
     bool found = false;
     while (it != takenRSide.end())
       {
	 if (*it == iTemp)
	   {
	     found = true;
	     break;
	   }
	 it++;
       }
     if (found)
       {
	 // This element is present
	   continue;
       }
     else 
       {
	 rIndex = iTemp;
	 break;
       }
	
   }
		// Stupid (for now) do the same stuff again
		
   for (int iTemp =0; iTemp < Ne_; iTemp++)
   {
     std::vector<int>::iterator it = takenLSide.begin();
     bool found = false;
     while (it != takenLSide.end())
       {
	 if (*it == iTemp)
	   {
	     found = true;
	     break;
	   }
	 it ++;
       }
     if (found)
       {
	 //	  This element is present
	   continue;
       }
     else 
       {
	 lIndex = iTemp;
	 break;
       }
     
   }
   return true;
}

/* Sets the aspect ratio (influences only SPwaveFct);
   if new_aToB<=0, nothing is changed. Returns 0 if change
   was successful and 1 if not. */
int yosState::setAspect(double new_aToB)
{
  if(new_aToB<=0) 
    {
      glLogger.info( "yosState::setAspect. Attempt to change aspect ratio to negative value or zero. Change denied");
      return 1;
    }
  aToB=new_aToB;
  return 0;
}

/*
  set fluxes of solenoids inside the torus
*/
void yosState::setSolenoidFluxes (double alpha1, double alpha2)
{
  alpha1_ = alpha1;
  alpha2_ = alpha2;
}

double yosState::getSolenoidFlux1() const {return alpha1_;}
double yosState::getSolenoidFlux2() const {return alpha2_;}



/* Single particle wavefunction of a state with momentum j. 
   Note that the WF is the same for spin up and spin down. */

std::complex<double> yosState::SPwaveFct(int j, 
				    double x, 
				    double y) const
{
  // x shift gives a y-dependent phase factor of exp
  std::complex<double> phaseFactor (1.0,0.0);
  if(y<0. || y>1.) 
    {
// move y into the unit cell
      glLogger.debug( "yosState::SPwaveFct called with y = (%f)  out of range [0;1]. y adjusted.", y);
if (y < 0.0)
{y = y+1.0;} else {y=y-1.0;}
      
    }
  if(x<0. || x>1.) 
    {
    	// move x into the unit cell
      glLogger.debug( "yosState::SPwaveFct called with x out of range [0;1]. x adjusted.");
            if (x < 0.0)
		{ x=x+1.0;}
else {x=x-1.0;}

	    phaseFactor = std::complex<double>(0.0,y * 2*M_PI* Nm_*aToB);
	    glLogger.debug("Phasefactor vor exp  is (%f), (%f)",
			   phaseFactor.real(), phaseFactor.imag()); 
	    phaseFactor = exp(phaseFactor);
	    glLogger.debug("Phasefactor is (%f), (%f)",phaseFactor.real(), phaseFactor.imag()); 
     }

  if (j < 0 || j >= Nm_) 
    {
     glLogger.info( "State number (j) is out of range [0..Nm-1]");
      return std::complex<double>(-1);
    }

  /*
    Declarations
  
    Names of variables are the same as in Chakraborty page 40
 
  */
  
  
  double Xj, ky, argImag, realExp, imagExp, prefact;

  /*
    cell is 1 x ab in units of the magnetic length i.e. l_[0] = 1
    without flexible  k-Sum for now, maybe later (i.e. k != 0)
  */
 
  realExp=0; imagExp=0;

  for (int k = -kSPsum ; k <= kSPsum; k++)
    {
      // j+1 for compatibility with fortran-code
      Xj = ( j + alpha2_/(2.0*M_PI) )/(1.*Nm_) + k;
      ky = 2*M_PI * ( j + Nm_*k );
      prefact = exp(-(Xj-x)*(Xj-x)*M_PI*Nm_*aToB);
      argImag = fmod(ky * y - alpha1_* (k + x), 2.0*M_PI);
      realExp = realExp + cos(argImag)*prefact;
      imagExp = imagExp + sin(argImag)*prefact;
    }

          /* To get really correct norm^{-1}, there should be
	     a factor 1/l_0, l_0=magnetic length. */

  std::complex <double> retVal =  std::complex <double>(realExp,imagExp) * sqrt(sqrt(2*aToB*Nm_))* phaseFactor;
  glLogger.debug("SPWavefct : Returning (%f), (%f)", retVal.real(), retVal.imag());
  return (retVal);
}

/* One particle wavefunction: derivative with respect to x of the j-th function at the place (x,y)
 * !!! this x is measured in units of cell with a, thus this is dphi/dx * a
 */
std::complex<double> yosState::SPwaveFctdx(int j,double x,double y) const
{
  double Xj, ky, alpha;

  std::complex<double> dphix (0.0, 0.0);
  /*
  for (int k = -kSPsum; k<= kSPsum; k++)
    {
      Xj = j/(1.*Nm_) + k;
      alpha = 2.0*M_PI * fmod( Xj * Nm_ * y, 1.0 ) - alpha1_*x - alpha2_*y;
      dphix +=  ( std::complex<double>(0.0, -alpha1_*x) + 2.0 * M_PI * Nm_ * aToB*(Xj-x) ) * 
	exp(-aToB*M_PI*Nm_*(Xj-x)*(Xj-x)) * std::complex<double>(cos(alpha), sin(alpha));
    }
  */


  for (int k = -kSPsum; k<= kSPsum; k++)
    {
      Xj = ( j + alpha2_/(2.0*M_PI) ) / (1.*Nm_) + k;
      ky = 2.0*M_PI * (j + Nm_*k); 
      alpha = fmod( ky * y - alpha1_* (k + x), 2.0*M_PI );
      dphix += ( 2.0*M_PI*Nm_*aToB * (Xj-x) - 
		 std::complex<double>(0.0, alpha1_) ) * exp(-aToB*M_PI*Nm_*(Xj-x)*(Xj-x)) * 
	std::complex<double>(cos(alpha), sin(alpha));
    }
  
  return sqrt(sqrt(2.0*aToB*Nm_)) * dphix;
}

/* One particle wavefunction: derivative with respect to y of the j-th function at the place (x,y)
 * !!! this y is measured in units of cell height b, thus this is dphi/dy * b* (a/b) = dphi/dy * a
 */
std::complex<double> yosState::SPwaveFctdy(int j,double x,double y) const
{
  double Xj, ky, alpha;

  std::complex<double> dphiy (0.0,0.0);
  /*
  for (int k= -kSPsum; k <= kSPsum; k++)
    {
      Xj = j/(1.*Nm_) + k;
      alpha = 2.0*M_PI * fmod( Xj * Nm_ * y, 1.0 ) - alpha1_*x - alpha2_*y;
      dphiy += (2.0*M_PI*Xj*Nm_ - alpha2_) * exp(-aToB*M_PI*Nm_*(Xj-x)*(Xj-x)) * 
	std::complex<double>(-sin(alpha), cos(alpha));
    }
  */

  for (int k= -kSPsum; k <= kSPsum; k++)
    {
      Xj = ( j + alpha2_/(2.0*M_PI) ) / (1.*Nm_) + k;
      ky = 2.0*M_PI * (j + Nm_* k);
      alpha = fmod( ky * y - alpha1_* (k + x), 2.0*M_PI );
      dphiy += ky * exp(-aToB*M_PI*Nm_*(Xj-x)*(Xj-x)) *	std::complex<double>(-sin(alpha), cos(alpha));
    }

  return aToB * sqrt(sqrt(2.0*aToB*Nm_)) * dphiy;
}



/*!
  \fn int * yosState::getjStates(void) const
  \return A pointer to the array of intger where the occ numbers are stored
*/


int * yosState::getjStates(void) const
{
  if (j_ != 0)
    {
      return j_;
    }
  else 
    {
    glLogger.info( "No states defined, returning 0 vector");
    return 0;
    }
}


/*!
  \fn int *yosState::getSpinStates(void) const
  \return A pointer to the array of intger where the occ numbers for the spin states are stored
*/

int *yosState::getSpinStates(void) const
{
  if (spinYes)
    {
      if (spin_ != 0)
	{
	  return spin_;
	}
      else 
	{
	  glLogger.info( "No spin values defined, return 0");
	  return 0;
	}
    }
  return 0;
}


int yosState::indInRange(int i) const
{
  if(i>=0 && i<Ne_) return i;
  else return 0;
}


/*! \fn int yosState::increaseNm(int increment) 
  \brief Increases Nm. Increment has to be positive.
 */
int yosState::increaseNm(int increment) 
{
  if(increment<0) return -1;
  Nm_ += increment;
  glLogger.debug("new Nm_:%d",Nm_);
  return 0;
}


void yosState::shiftRight_State(int shift)
{
  for(int i=0;i<Ne_;i++)
    j_[i]=((j_[i] + shift) % Nm_);
}

/*! \fn int yosState::compareAndFindDiffs(const yosState& rhs,
				  cmp_yosStates_oneDifference &diff)
    \brief Find differences between this state and rhs. Provided
    there are max. MAX_DIFF different 1-part. states describe the difference in
    diff.
 */
int yosState::compareAndFindDiffs(const yosState& rhs,
				  cmp_yosStates_differences &diff) const
{
  if(Nm_ != rhs.getNm() || Ne_ != rhs.getNe()) 
  {
    glLogger.warning("yosState::compareAndFindDiffs");
    glLogger.warning("comparing two states of unlike Ne/Nm type");
    return Ne_ + rhs.getNe();
  }

  /*
  printf("\n");
  print();
  rhs.print();
  */

  static const int MAX_DIFF = cmp_yosStates_differences::MAX_DIFF;

// Go through states in lhs and see if they are also present in rhs
  int i_r=0;
  int i_r_orig;
  int diffs_l=0;
  for(int i_l=0;i_l<Ne_;i_l++,i_r++)
  {
    if(st_nr(i_l) == rhs.st_nr(i_r) && i_r<Ne_) continue;
    i_r_orig = i_r;
    do i_r++;
    while(i_r<Ne_ && st_nr(i_l) != rhs.st_nr(i_r));
    if(i_r >= Ne_) // j_[i_l],spin_[i_l] is not in rhs
    {
      diff.pos_diff_l [diffs_l % MAX_DIFF]  = i_l;
      diff.j_diff_l   [diffs_l % MAX_DIFF] = j_[i_l];
      diff.spin_diff_l[diffs_l % MAX_DIFF] = spin_[i_l];
      i_r = i_r_orig - 1;
      diffs_l++;
    }
  }

// Go through states in rhs and see if they are also present in lhs
  int i_l=0;
  int i_l_orig;
  int diffs_r=0;
  for(int i_r=0;i_r<Ne_;i_l++,i_r++)
  {
    if(st_nr(i_l) == rhs.st_nr(i_r)) continue;
    i_l_orig = i_l;
    do i_l++;
    while(i_l<Ne_ && st_nr(i_l) != rhs.st_nr(i_r));
    if(i_l >= Ne_) // j_[i_l],spin_[i_l] is not in rhs
    {
      diff.pos_diff_r [diffs_r % MAX_DIFF] = i_r;
      diff.j_diff_r   [diffs_r % MAX_DIFF] = rhs.j(i_r);
      diff.spin_diff_r[diffs_r % MAX_DIFF] = rhs.spin(i_r);
      i_l = i_l_orig - 1;
      diffs_r++;
    }
  }
  
  /*
  printf("(%d,%d:",diffs_l,diffs_r);
  for(int i=0;i<diffs_l;i++)
    printf("%d,",diff.pos_diff_l[i]);
  printf(";");
  for(int i=0;i<diffs_r;i++)
    printf("%d,",diff.pos_diff_r[i]);
  printf(")");
  */

  if(diffs_r != diffs_l) 
  {
    glLogger.error("yosState::compareAndFindDiffs (twoDifferences)");
    glLogger.error("comparison algorithm wrong (diffs_l!=diffs_r)");
    exit(1);
  }
  return diffs_r;
}

void yosState::sPrint_num(char str[],int len) const
{
  int pos=0;

  if(len<2)
    {
      glLogger.error("yosState::sPrint_num");
      glLogger.error("Too short string passed.");
      exit(1);
    }
  sprintf(str,"|");
  pos++;
  for(int i=0;i<Ne_;i++)
    {
      if(pos+5>len)
	{
	  glLogger.error("yosState::sPrint_num");
	  glLogger.error("Too short string passed.");
	  exit(1);
	}

      sprintf(str+pos,"%d",j_[i]);
      pos ++;
      if(j_[i]>9) pos ++;
      if(j_[i]>99) pos ++;
      
      if(spinYes)
	{
	  if(spin_[i]) sprintf(str+pos,"+");
	  else
	    sprintf(str+pos,"-");
	  pos++;
	}
    }
  sprintf(str+pos,">");
}

void yosState::sPrint_occ(char str[],int len) const
{
  if(len<Nm_+3)
    {
      glLogger.error("yosState::sPrint_occ");
      glLogger.error("Too short string passed.");
      exit(1);
    }  

  for(int i=0;i<Nm_+2;i++) str[i]=32;  // 32 = space
  str[0] = 124;  // 124 = |
  // 23 = up, 24 = down
  for(int i=0;i<Ne_;i++)
    {
      if(spinYes)
	{
	  if(str[j_[i]+1]==32)
	    if(spin_[i])
	      str[j_[i]+1]=94;  // 94 = ^ (up) 
	    else 
	      str[j_[i]+1]=118;  // 118 = v (down)
	  else
	    str[j_[i]+1]=73;    // 73 = I (i.e. up and down)
	}
      else
	str[j_[i]+1]=120;  // 120 = x
    }
  str[Nm_+1]=62;
}




/*! \fn int yosState::st_nr(int i)
  \brief A way to characterize the i-th one-particle state
  of this state by one number uniquely (2*j+spin). This allows 
  for a simpler (faster) way to find out whether two (j,spin) 
  one-particle states are identical or not.
 */
int yosState::st_nr(int i) const
{
  return 2*j_[i % Nm_]+spin_[i % Nm_];
}


