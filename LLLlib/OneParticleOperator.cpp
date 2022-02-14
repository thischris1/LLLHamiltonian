#include <OneParticleOperator.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>

#include <fstream>
#include <iostream>


  ////////////////////
 // static members //
////////////////////

// needed for sharing landau marices among several one-particle-ops
int OneParticleOperator::m_iTypeOfLandauMatrix = 0;
bool OneParticleOperator::m_bLandauEvaluated = false;
eigSt OneParticleOperator::m_stateLandauEvaluated;
const std::complex<double> *OneParticleOperator::m_pGlobalLandauMatrix = NULL;


  /////////////
 // members //
/////////////

OneParticleOperator::OneParticleOperator(const Basis &newBasis, 
					 const eigSt &state, 
					 int new_type,
					 bool bNeedsWF,
					 bool bDoNotComputeLandauMatrix) : 
  Operator( (Basis*) &newBasis, new_type)
{
  m_bNeedsWF = bNeedsWF;
  
  m_stateToCompute = state;
  m_bStateValid = true;

  landau_matrix_allocated = false;


  //  std::cout << "constructor OneParticleOperator::OneParticleOperator(const Basis &newBasis, 
  //				 const eigSt &state, 
  //				 int new_type) \n";

  endYosBasis=(yosBasis *)endBasis;
  Nm = endYosBasis->getNm();
  Ne = endYosBasis->getNe();

  /* This operator can handle only (linear combinations of)
     yosStates, i.e. endBasis has to be of type yosBasis. */
  if(endBasis->getBasisType()!=YOS_BASIS) {
  	throw (CxBadValueError(__FILE__,__LINE__, "DensOperator may be used only for bases: yosBasis or linear combination of such basis. Other types must be implemented."));
  
  }

 
  if (bNeedsWF)
    endYosBasis->initWFCache();

  /* Compute the Landau matrix (LM); 
     omitting this makes only sense if the ctor of a derived
     class calls setupLandauMatrixForState().
  */
  if(!bDoNotComputeLandauMatrix) setupLandauMatrixForState ();  

}

OneParticleOperator::OneParticleOperator( const Basis &newBasis, 
					  const OneParticleOperator &ancestor, 
					  int new_type,
					  bool bNeedsWF) : 
  Operator( (Basis*) &newBasis, new_type),
  st1(0),
  st2(0),
  Nm(),
  Ne(),
  x(0.0),
  y(0.0)
{
  assert (false); // not yet tested...

  landau_matrix_allocated = false;

  m_bStateValid = false;
  m_bNeedsWF = bNeedsWF;

  endYosBasis=(yosBasis *)endBasis;
  Nm = endYosBasis->getNm();
  Ne = endYosBasis->getNe();

  /* This operator can handle only (linear combinations of)
     yosStates, i.e. endBasis has to be of type yosBasis. */
  if(endBasis->getBasisType()!=YOS_BASIS){
  	throw CxBadValueError(__FILE__,__LINE__);
  }

  if (bNeedsWF) 
    endYosBasis->initWFCache();
 

  m_pCmplxLandauMatrix = allocK< std::complex<double> >(Nm*Nm);
  landau_matrix_allocated = true;

  memcpy ( m_pCmplxLandauMatrix, ancestor.getLandauMatrix(), Nm*Nm*sizeof(std::complex<double>) );

}

OneParticleOperator::~OneParticleOperator()
{
  if(landau_matrix_allocated)
    {
      delete [] m_pCmplxLandauMatrix;
      landau_matrix_allocated=false;
    }
}

int OneParticleOperator::preprocessMatEls(int newBasis_matrix_type)
{
  int dim,dim1;
  std::complex<double> xtmp;

  if(landau_matrix_allocated) return 0;  // Already generated.

  switch(my_basis->getBasisType()){
  case LINEAR_COMBINATION:
    {
      assert (m_bStateValid); 
      assert(false); /* This is an implementation done at time when there was only
			densOperator (and not OneParticleOperator, spinDensOperator,
			etc.) So as it would work, operator of a correct type has to
			be instantiated below (i.e. if you call this method
			from spinDensOperator, op has actually to be of 
			spinDensOperator type). */
      
    if(newBasis_matrix_type!=0) 
    {
    	throw CxBadValueError(__FILE__,__LINE__);
    }
    OneParticleOperator op(*my_basis->getMyBasis(), m_stateToCompute, newBasis_matrix_type, false,false);
    op.preprocessMatEls(newBasis_matrix_type);
    dim=my_basis->dimension();
    dim1=my_basis->getMyBasis()->dimension();    

    m_pCmplxLandauMatrix = allocK< std::complex<double> >(dim*dim*Nm*Nm);

    landau_matrix_allocated=1;
    // Compute the matEl between i-th and j-th vector of the basis
    for(int i=0;i<dim;i++)     // (LC) basis vectors
      for(int j=i;j<dim;j++)
	{
	  for(int m=0;m<Nm;m++)    // Elements of the Landau matrix
	    for(int n=0;n<Nm;n++)
	      {
		xtmp = std::complex<double>(0.0,0.0);
		for(int i1=0;i1<dim1;i1++)  // End-basis vectors
		  for(int j1=0;j1<dim1;j1++)
		    xtmp = xtmp + op.allLandauMatricesElement(i1,j1,m,n)*my_basis->getCoef(i,i1)*my_basis->getCoef(j,j1);
		m_pCmplxLandauMatrix[i*dim*Nm*Nm+j*Nm*Nm+m*Nm+n] = xtmp;
	      }
	}
    return 1;
     
    }
  case YOS_BASIS:  // No preprocessing needed: just ask for matEls
    setupLandauMatrix();
    return 1;
  default:
   throw CxErrors(__FILE__,__LINE__);
  }
  return 0;
}

std::complex<double> OneParticleOperator::allLandauMatricesElement( int i_st1, int i_st2, int i, int j )
{
  int dim;
  if(landau_matrix_allocated==0) setupLandauMatrix();
  if(landau_matrix_allocated==0) {
  	throw CxErrors("Got no matrix", __FILE__,__LINE__);
  }
  dim=my_basis->dimension();
  if(i_st1<0 || i_st1>=dim) {
  	throw CxErrors(__FILE__,__LINE__);
  }
  if(i_st2<0 || i_st2>=dim) {
  	throw CxErrors(__FILE__,__LINE__);
  }
  if(i<0 || i>=Nm) {
 
  	throw CxErrors(__FILE__,__LINE__);

  	}
  if(j<0 || j>=Nm) {
  	throw CxErrors(__FILE__,__LINE__);
  	}
  return m_pCmplxLandauMatrix[i_st1*dim*Nm*Nm+i_st2*Nm*Nm+i*Nm+j];
}




/* This function computes the matrix elements (i.e. a Landau matrix for each 
   pair of basis states) in yosBasis.
   "AllLandauMatrix"
*/


int OneParticleOperator::setupLandauMatrix(void)
{
  int dim = endBasis->dimension();
  if(endBasis!=my_basis) {
  	
  	throw CxErrors(__FILE__, __LINE__);
  }
  m_pCmplxLandauMatrix = allocK< std::complex<double> >(dim*dim*Nm*Nm);
  landau_matrix_allocated=1;
  for(int i1=0;i1<dim;i1++)    // Run through all the basis states
    for(int i2=0;i2<dim;i2++)
      {
	for(int m=0;m<Nm;m++)
	  for(int n=0;n<Nm;n++)
	    m_pCmplxLandauMatrix[i1*dim*Nm*Nm+i2*Nm*Nm+m*Nm+n] = std::complex<double>(0.0,0.0);
	  addToLandauMatrix(i1, i2, m_pCmplxLandauMatrix+i1*dim*Nm*Nm+i2*Nm*Nm, 1.0);
      }
  return 0;
}


// precondition: object must be constructed through constructor
// that defines the state, not the one that inherits a landaumatrix
int OneParticleOperator::setupLandauMatrixForState()
{
  assert (m_bStateValid);

  // first check, if landau-matrix was already calculated by another instance
  // of a OneParticleOpertor derived class
  if ( m_bLandauEvaluated && 
       (m_stateToCompute == m_stateLandauEvaluated) &&
       m_iTypeOfLandauMatrix == whatTypeOfLandauMatrixDoINeed() )
    {
      // draw a copy of already (by another one-particle-op) evaluated landau matrix
      glLogger.info ("Found valid landau matrix of another operator (Landau matrix type %d; my own type %d.", m_iTypeOfLandauMatrix, whatTypeOfLandauMatrixDoINeed());

      m_pCmplxLandauMatrix = allocK< std::complex<double> >(Nm*Nm);
      landau_matrix_allocated = true;
      memcpy ( (void*) m_pCmplxLandauMatrix, m_pGlobalLandauMatrix, Nm*Nm*sizeof(std::complex<double>) );
    }
  else
    {
      // has to be calculated

      /* Matrix elements (actually, each element is a Landau
	 matrix) of the density operator will be stored
	 as FULL_REAL. */
      //int new_matrix_type=0;   

      int dim = endBasis->dimension();

      if( !landau_matrix_allocated )
	{
	  m_pCmplxLandauMatrix = allocK< std::complex<double> >(Nm*Nm);
	  landau_matrix_allocated=true;

	  //preprocessMatEls(new_matrix_type);
      
	  // Clean up the Landau matrix
	  for(int m=0;m<Nm;m++)
	    for(int n=0;n<Nm;n++)
	      m_pCmplxLandauMatrix[m*Nm+n] = std::complex<double>(0.0,0.0);
      
	  // Set it up (correspondingly to coefficients in state)
	  glLogger.info("Setting up Landau matrix (only lower triangle, m(i,j), j<i) (LM type %d)", whatTypeOfLandauMatrixDoINeed());
	  for(int i1=0; i1<dim; i1++)  // End-basis vectors
	    {

	      for(int i2=0; i2<=i1; i2++)
	        addToLandauMatrix(i1, i2, m_pCmplxLandauMatrix, conj(m_stateToCompute[i1]) * m_stateToCompute[i2]);

	      //for(int i2=0; i2<dim; i2++)
	      //addToLandauMatrix(i1, i2, m_pCmplxLandauMatrix, conj (m_stateToCompute[i1]) * m_stateToCompute[i2]);

	      // progress-logging
	      if(!(i1%FREQ_REPORT_OPOP))
		{
		  if(!i1) Npercent(-1.,stdout);
		  Npercent(1.*i1/dim,stdout);
		}
	    }
	  glLogger.info (" done.");
     
	  // Print the Landau matrix. It should be symmetric.

	  glLogger.debug ("Landaumatrix = ");
	  for(int i=0;i<Nm;i++)
	    for(int j=0;j<Nm;j++)
	      {
		std::complex<double> z = allLandauMatricesElement(0,0,i,j);
		glLogger.debug ( "L(%d,%d)=(%f,%f)", i, j, z.real(), z.imag() );
	      }


	}

      // remember for which state we computed the landau-matrix
      // and initialize globally accessable const pointer
      if(!m_bLandauEvaluated)
	{
	  m_pGlobalLandauMatrix = m_pCmplxLandauMatrix;
	  m_bLandauEvaluated = true;
	  m_iTypeOfLandauMatrix = whatTypeOfLandauMatrixDoINeed();
	  m_stateLandauEvaluated = m_stateToCompute;         
	}
    }

  return 0;
}

/* Computes the Landau matrix between i1-th and i2-th state of the 
   endYosBasis and adds the result (multiplied by fac) to array 
   pointed to by 'where'. Returns zero if nothing was changed ('matrix
   element' is zero) and one otherwise. 
   
   Note: i1 should be <= i2; (i2,i1) gives the same contribution to 
   the Landau matrix as (i1,i2) but places it into (j2diff,j1diff)-element
   instead of (j1diff,j2diff)-element. See below.
*/ 
int OneParticleOperator::addToLandauMatrix(int i1,int i2, std::complex<double> *where, std::complex<double> fac)
{
  int j1,j2,k1,k2,j1diff,j2diff,j_i;
  std::complex<double> newMatEl;
  
  st1=(*endYosBasis)[i1];
  st2=(*endYosBasis)[i2];
  //printf("Comparing: ");st1->print();st2->print();printf(": ");
  // Compare the two states
  j1=0;j2=0;
  if(i1==i2) j1=Ne; // The states are definitely the same
  while(s1_nr(j1)==s2_nr(j2) && j1<Ne) {j1++;j2++;}
  if(j1<Ne)
    {
      if(j1==Ne-1) // and thus j2=Ne-1;   | ..345>, | ..347>
	{
	  k1=j1;k2=j2;j1diff=s1_nr(j1)/2;j2diff=s2_nr(j2)/2;
	}
      else
	if(s1_nr(j1+1)==s2_nr(j2)) // |..3456..>, |..356..>
	  {
	    k1=j1;j1diff=s1_nr(j1)/2;
	    j1++;
	    while(s1_nr(j1)==s2_nr(j2) && j1<Ne) {j1++;j2++;}
	    j2diff=s2_nr(j2)/2;k2=j2;
	    j2++;
	    while(s1_nr(j1)==s2_nr(j2) && j1<Ne) {j1++;j2++;}	    
	    if(j1<Ne) return 0; /* The states differ in more than one QN,
				 => no contribution to LandauMatrix.*/
	  }
	else
	  if(s1_nr(j1)==s2_nr(j2+1)) // |..356..>, |..3456..>
	    {
	      k2=j2;j2diff=s2_nr(j2)/2;
	      j2++;
	      while(s1_nr(j1)==s2_nr(j2) && j2<Ne) {j1++;j2++;}
	      j1diff=s1_nr(j1)/2;k1=j1;
	      j1++;
	      while(s1_nr(j1)==s2_nr(j2) && j1<Ne) {j1++;j2++;}	    
	      if(j2<Ne) return 0; /* The states differ in more than one QN,
				     => no contribution to LandauMatrix.*/
	    }  
	  else                       // |..3467..>, |..3567..>
	    {  
	      k1=j1;k2=j2;j1diff=s1_nr(j1)/2;j2diff=s2_nr(j2)/2;
	      j1++;j2++;
	      while(s1_nr(j1)==s2_nr(j2) && j1<Ne) {j1++;j2++;}
	      if(j1<Ne) return 0; /* The states differ in more than one QN,
				     => no contribution to LandauMatrix.*/
	    }
      /* The two multipart(MP)-states also differ in one 1-part state;
	 in the 1st MPstate on the place k1 there's a 1-part state
	 (with j=j1diff) which is not in the 2nd MPstate.
	 in the 2nd MPstate on the place k2 there's a 1-part state
	 (with j=j2diff) which is not in the 1st MPstate. */
      //printf("left at %d is %d, right at  %d is %d.\n",k1,j1diff,k2,j2diff);fflush(stdout);

      newMatEl = fac*getOneElementOfLandauMatrix(k1,k2,i1,i2);

      //printf("/=%f\n",newMatEl/fac);
      if(((k1+k2)%2)==1) newMatEl=-newMatEl;
      /* Asuming that addToLandauMatrix() is called only for 
	 (i1,i2) and not for (i2,i1). */
      where [j1diff*Nm+j2diff] += newMatEl;
      where [j2diff*Nm+j1diff] += conj(newMatEl);
      return 1;
    }
  /* The two multipart-states are the same */
  if(i1!=i2) fac=fac*2.; /* For similar reasons as above 
			   (j1diff,j2diff) vs. (j2diff,j1diff) */
  for(int i=0; i<Ne; i++)
    {
      j_i=st1->j(i);
      newMatEl = fac*getOneElementOfLandauMatrix(i,i,i1,i2);
      where[j_i*Nm+j_i] += newMatEl;
    }

  return 1;
}


/* This method is expected to be overloaded in derived classes. 
   It is used in setupLandauMatrix() only (and see  */
std::complex<double> OneParticleOperator::getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2)
{
/* The two states differ by the k1-th one-el. wavefunction of st1 and 
the k2-th wavefunction of st2; is these two one-el. wavefunctions have 
different spin then <..|\delta(r-r')|..>=0 because \delta(..) has no spin part. 

i1 and i2 are the indices of the two yosStates between which the matEl is 
computed. They are not needed exactly here but are useful in derived classes.

This condition should be actually tested only if st1!=st2. For st1=st2 this 
routine should thus be called in a manner which ensures that the following if()
condition is _not_ fulfilled (i.e. return 1.0).
*/
  if((s1_nr(k1) % 2) != (s2_nr(k2) % 2)) 
    return std::complex<double>(0.0,0.0);
  else
    return std::complex<double>(1.0,0.0);
}

int OneParticleOperator::whatTypeOfLandauMatrixDoINeed()
{ return 1;}


int OneParticleOperator::setXY(double new_x, double new_y)
{
  //std::cout << "setXY\n";

  int retval=0;
  if(new_x<0. || new_x>1.) 
    {
      std::cerr << "OneParticleOperator::setXY called with x out of range [0;1].\n";
      new_x=(new_x<0 ? 0. : 1.);
      retval=1;
    }
  if(new_y<0. || new_y>1.)
    {
      std::cerr << "OneParticleOperator::setXY called with y out of range [0;1].\n";
      new_y=(new_y<0 ? 0. : 1.);
      retval=1;
    }
  
  if (m_bNeedsWF)
    endYosBasis->setWFCacheForXY(new_x, new_y);

  x = new_x;
  y = new_y;

  //std::cout << "setXY done\n";      
  return retval;
}

int OneParticleOperator::matEl_nonzero(int i_st1,int i_st2, std::complex<double> *a)
{
  std::complex<double> Oij;

  if(i_st1>=my_basis->dimension() || i_st1<0) 
  {
  	
  	throw CxErrors(__FILE__, __LINE__);
  }
  if(i_st2>=my_basis->dimension() || i_st2<0) {
  	throw CxErrors(__FILE__, __LINE__);
  }

  if(!landau_matrix_allocated) setupLandauMatrixForState ();  
  //preprocessMatEls(new_matrix_type);
  *a=std::complex<double>(0,0);

  for(int i=0; i<Nm; i++)
    {
      (*a) += allLandauMatricesElement(i_st1,i_st2,i,i) * getOneParticleEl (i, i);
      for(int j=0;j<i;j++)
	{
	  Oij = getOneParticleEl (i, j);	  
	  *a += allLandauMatricesElement(i_st1,i_st2,i,j) * Oij + allLandauMatricesElement(i_st1,i_st2,j,i) * conj(Oij);
	}
    }

  return 1;   // Yes, this matEl is non-zero... (actually, this is not that sure)
}

double OneParticleOperator::evalAt (double x, double y)
{
  std::complex<double> val (0.0,0.0), Oij;

  setXY (x, y);

  // assuming symmetric Landau-Matrix

  if(!landau_matrix_allocated) setupLandauMatrixForState ();  

  for(int i=0; i<Nm; i++)
    {
      val += allLandauMatricesElement(0,0,i,i) * getOneParticleEl (i, i);
      for(int j=0; j<i; j++)
	{
	  Oij = getOneParticleEl (i, j);
	  //val += allLandauMatricesElement(0,0,i,j) * ( Oij + conj(Oij) );
	  val += allLandauMatricesElement(0,0,i,j) * Oij + allLandauMatricesElement(0,0,j,i) * conj(Oij);
	}
    }

  return real(val);
}



inline int OneParticleOperator::s1_nr(int i) const
{return (int) st1->j(i)*2+st1->spin(i);}
inline int OneParticleOperator::s2_nr(int i) const
{return (int) st2->j(i)*2+st2->spin(i);} 





















