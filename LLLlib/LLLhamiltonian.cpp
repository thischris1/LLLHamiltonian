/*! 
  \file LLLhamiltonian.cpp
 */

/*!
  \class LLLhamiltonian : public Operator
  \brief Implementation of Hamiltonian of the Lowest Landau Level, i.e.
  system of Ne electrons in a rectangle with periodic boundary conditions which
  sit in the LLL and interact by Coulomb interaction.

  This is a homogeneous system, if inhomogeneities are desired, derive a class
  from LLLhamiltonian and rewrite matEl_perturb_nonzero(i_st1,i_st2,a) 
  (this gives the matrix element (a) of the inhomogeneity between i_st1-th and
  i-st2-th state of the basis). See for example LLLhamInplaneMagImp.

  Between constructing and accessing the matrix elements via matEl_nonzero you  should call
  either computeCoulomb() or computeHardcore() (which computes the
  2-particle matrix elements).

*/


#include<LLLhamiltonian.hpp>
#include<th1_natureConst.hpp>
#include <assert.h>
#include<math.h>
#include<cstdlib>

#include <utils/logger.hpp>
#include <ERRORS.h>

#ifdef NAG
#else
#include <gsl/gsl_sf_gamma.h>
#endif
/*
  static members
*/
const double LLLhamiltonian::CUTOFF_MIN_Q2=1.e-6;   // in computeCoulomb: cutoffs for
const double LLLhamiltonian::CUTOFF_MAX_Q2=50.e0;   // q=(q_x,q_y) (in MIN<|q|^2<MAX)

/*! 
  \fn LLLhamiltonian::LLLhamiltonian(yosBasis* new_basis,int new_type,double bli_new, double a_new, double b_new)

  \param new_basis The basis with respect to which 
                   we will be calculating the matrix elements.
  \param new_type  How to store the matrix elements (0=FULL_REAL,
                        1=SPARSE_REAL, 2=FULL_COMPLEX); see also Operator.cpp
  \param bli_new, a_new, b_new Parameters of the system ("physics"; bli =
                   finite width parameter, a,b size of the system; note
		   that only a:b=aspect ratio matters, the real size of the 
		   system is given by Ne, Nm which come from the basis).
*/

LLLhamiltonian::LLLhamiltonian(yosBasis* new_basis,int new_type,double bli_new, 
								double a_new, double b_new) 
: Operator(new_basis,new_type)
,LLL_basis(new_basis),
Ne(),
Nm(new_basis->getNm()),
bli(bli_new),
a(a_new),
b(b_new),
coulombME(0)
{
  double scal;   // see below
  if (!LLL_basis)
  {
  	throw CxNullPointerError(__FILE__,__LINE__);
  }
  Ne= new_basis->getNe();
  Nm= new_basis->getNm();
  
  //basis_type=YOS_BASIS;
  
  switch(new_basis->getSpinYes())
  {	
  	case 0: 
  		spinPolarized=1;
  				break;
	  default: 	
	  	spinPolarized=0;
  }

  scal=sqrt(a*b/2/M_PI/Nm);   // reads ab/(2*pi*Nm)=1. Scaling done here
  a=a/scal;
  b=b/scal;   // keeps a/b and makes the mentioned cond. be fulfilled.
  coulombME_computed=0;

  coulombME=new double*[Nm];
  
  for(int i=0;i<Nm;i++) 
    {      
      coulombME[i]=new double[Nm];
          }

  //if(computeCoulomb()) err_msg(9);  Moved to the public part...
}

LLLhamiltonian::~LLLhamiltonian()
{
  if(coulombME)
    {
      for(int i=0;i<Nm;i++) delete[] coulombME[i];  
      delete[] coulombME;
      coulombME = 0;
    }
}


/*! \fn double LLLhamiltonian::computeMadelung() 

  \brief The Madelung (i.e. 1-particle) part of the Hamiltonian ('Wigner
  crystal energy' - interaction of electron with its periodic images). 
  This is actually just an energy shift.
  \returns The energy shift due to interaction of electrons with their 
  periodic images.

*/

double LLLhamiltonian::computeMadelung()
{
  double s=0;    // Here should be the sum with phi_{-1/2}, Chakr. p44
  const static double pi=M_PI;
  double al=0.5;
double minarg = 0.0;
double x = 0.0;
double arg = 0.0;
  
	int 
  	l1 = 0,
  	l2 = 0,
  	ifail = 0;
  const static double MAX_GAMMAFCT_ARG=20.0; // for x=12, G(1/2,x) is about 10^-6
  l1=0;
  while(minarg<MAX_GAMMAFCT_ARG)
    {
      l2=0;
      if(l1==0) l2=1;
      while((arg=pi*(a/b*l1*l1+b/a*l2*l2))<MAX_GAMMAFCT_ARG || l2<2)
	{
	  /* see the NAG routine documentation; 
	     \phi_{-1/2}(z)=\sqrt(\pi/z)*Q(1/2,z) */
	  #ifdef NAG

		  s14baf_(&al,&arg,&PREC_GAMMAFCT,&foo,&x,&ifail);
	  #else
	  glLogger.debug("Call Gsl gamma");
	  x = gsl_sf_gamma_inc_Q (al, arg);
	   glLogger.debug("Retrieved result (%f)", x);
	  #endif
	  if(ifail) 
	  {
	  	throw (CxBadValueError(__FILE__,__LINE__," ifail after minmax function"));
	  	
	  }
	  x=x*sqrt(pi/arg);
	  glLogger.debug("x in computeMadelung is   (%f)", x);
	  if((l1==0)||(l2==0)) s+=x/2.;  // count (0,l), (l,0) twice
	  else s+=x;       // and the rest four times (+-l1,+-l2)
	  if(l2==0) minarg=arg;
	  l2++;
	}
      l1++;
    }
  s=4*s;
  s = -(2.-s)/sqrt(2*pi*Nm);
  glLogger.info(" Madelung Energie %f", s);
  return s;   // Chak., section 'Rectangular geom.'
}

/*! \fn int LLLhamiltonian::computeCoulomb()
  \brief Computes the 2-particle matrix elements of Coulomb interaction and 
  stores them to coulombME[][]. At the end computeMadelung() is called.

  The many-particle matrix elements are then computed from coulombME by 
  matEl_nonzero.
*/
int LLLhamiltonian::computeCoulomb()
{
  double abl,cf,q,q2,b1,b2;
  int km1,km2;
  const double pi=M_PI;

///////////////////////// Is really t=j1-j4? (see Chak., spin paper PRB 30)

  glLogger.info(" - Computing the two-particle operator matrix elements");

  b1=sqrt(2*pi*Nm*a/b);
  b2=sqrt(2*pi*Nm*b/a);
  km1=(int)(ceil(b1))+1;   // maximal k1,k2; these correspond to 
  km2=(int)(ceil(b2))+1;   // q=(k1,k2) with q^2>=2*pi*l_0^2,
  // which contributes negligibly ( prop. to exp(-q^2 l_0^2) )
  
  // Particular matrix elements A_{j1,j2,j3,j4}; ja=j1-j4, jb=j1-j3.
  // A_{j1,j2,j3,j4}=A_{ja,jb}; note that A_{ja,jb}=A_{ja,-jb},
  // A_{ja,jb}=A_{-ja,jb} (only if the elements remain real)
  // !!!! NOTE: ja,jb=1,..,Nm in fqper_sparV and ja,jb=0,..,Nm-1 here !!!!

  for(int ja=0;ja<Nm;ja++)    
    for(int jb=0;jb<Nm;jb++)  
      {
	coulombME[ja][jb] = 0.0;
        for(int k1=-km1;k1<km1+1;k1++)
	  for(int k2=-km2;k2<km2+1;k2++)
	    {
	      if((ja-k2)%Nm) continue;  // proceed only if (ja-k2 mod Nm)=0
	      q2=pow(2*pi/b1*k1,2)+pow(2*pi/b2*k2,2); // q_x^2+q_y^2
	      if(q2<CUTOFF_MIN_Q2 || q2>CUTOFF_MAX_Q2) continue;
	      q=sqrt(q2);
	      // See Chakraborty, p.66, (5.64)
	      abl=(8.e0+9.e0*q*bli+3.e0*pow(q*bli,2))/pow(1.e0+q*bli,3)/8.e0;
	      // See Chakraborty, p.41, (5.6); CAUTION: real part of it
	      cf=exp(-q2*5.e-1)*cos(2*pi*k1*jb/Nm);
	      coulombME[ja][jb]+=cf/q*abl/2.e0/Nm;
	    }
      }
  madelung=computeMadelung();           

  glLogger.info("done.");
  // Print the matrix elements to the log file.
  
  glLogger.debug("'Coulomb' (two-particle operator) matrix elements, A_{j1,j2,j3,j4}=A_{ja,jb}, ja=j1-j4, jb=j1-j3.\n\n");
  for(int i=0;i<Nm;i++)
    {
      for(int j=0;j<Nm;j++)
	{
	  glLogger.debug("(%d), (%d) %8.3g ",i, j, coulombME[i][j]);
	}
    }
      
  coulombME_computed=1;
  return 0;
}

/*! \fn int LLLhamiltonian::computeHardcore()
  \brief Similar like computeCoulomb() but asumes hardcore interaction instead
  of Coulomb interaction. Use either computeHardcore() or computeCoulomb(), 
  not both.
*/
int LLLhamiltonian::computeHardcore()
{
  double q2,b1,b2, sumK2;
  int km1,km2;
  const double pi=M_PI;
  double k0x2, k0y2; 


  glLogger.error(" - Computing the two-particle operator matrix elements of hardcore-Interaction");
  
  b1=sqrt(2*pi*Nm*a/b);
  b2=sqrt(2*pi*Nm*b/a);
  k0x2 = pow (2*pi/b1, 2);
  k0y2 = pow (2*pi/b2, 2);
  

  km1=(int)(4*(ceil(b1)+1));   // maximal k1,k2; these correspond to 
  km2=(int)(4*(ceil(b2)+1));   // q=(k1,k2) with q^2>=2*pi*l_0^2,
 
  /*
    which contributes negligibly ( prop. to exp(-q^2 l_0^2) )
  
   Particular matrix elements A_{j1,j2,j3,j4}; ja=j1-j4, jb=j1-j3.
   A_{j1,j2,j3,j4}=A_{ja,jb}; note that A_{ja,jb}=A_{ja,-jb},
   A_{ja,jb}=A_{-ja,jb} (only if the elements remain real)
   !!!! NOTE: ja,jb=1,..,Nm in fqper_sparV and ja,jb=0,..,Nm-1 here !!!!

  */

  for(int ja=0; ja<Nm; ja++)    
    for(int jb=0; jb<Nm; jb++)  
      {
	coulombME[ja][jb] = 0.0;

        for(int k1=-km1; k1<=km1; k1++)
	  {
	    sumK2 = 0.0;
	    for(int k2=-km2; k2<=km2; k2++)
	      {
		if( (ja-k2) % Nm ) continue;  // proceed only if (ja-k2 mod Nm)=0

		q2 = k0x2 * k1*k1 + k0y2 * k2*k2; // q_x^2+q_y^2

		//if( q2 > CUTOFF_MAX_Q2 ) continue;

		sumK2 += exp(-q2*5.e-1) * q2;
	      }
	    coulombME[ja][jb] -= sumK2*cos(2*pi*k1*jb/Nm)/(2.e0*Nm);
	  }
      }
/*
 * No Madelung energy for Short range energy 
 * But does however only contribute a constant value to the energy (0.20087 for 5 electrons)
 glLogger.error("Here");
  madelung=computeMadelung();           
  */
  glLogger.debug(", Madelung=%f done",madelung);


  // Print the matrix elements to the log file.
  if(logfile!=NULL)
    {
      fprintf(logfile,"\n'Hardcore' (two-particle operator) matrix elements, A_{j1,j2,j3,j4}=A_{ja,jb}, ja=j1-j4, jb=j1-j3.\n\n");
      for(int i=0;i<Nm;i++)
	{
	  for(int j=0;j<Nm;j++)
	    fprintf(logfile,"%8.3g ",coulombME[i][j]);
	  fprintf(logfile,"\n");
	}
    }

  coulombME_computed=1;

  return 0;
}

/*!
 *
 *
 */
int LLLhamiltonian::computeNoInteraction()
{
	for(int ja=0; ja<Nm; ja++) {
	    for(int jb=0; jb<Nm; jb++)
	      {
		coulombME[ja][jb] = 0.0;
	      }
	}

	coulombME_computed=1;

	  return 0;

}
/*! \fn int LLLhamiltonian::matEl_perturb_nonzero(int i_st1,int i_st2,std::complex<double> *a)

\brief Evaluation of matrix element of the perturbation should be implemented here
returns zero if the ME is zero and one otherwise. The value of ME is in a.
*/
int LLLhamiltonian::matEl_perturb_nonzero(int ,int ,std::complex<double> *a)
{

  (*a)=std::complex<double>(0.0,0.0);
  return 0;
}

/*! \fn std::complex<double> LLLhamiltonian::matEl_perturb(int i_st1,int i_st2)
    
\brief Returns just the value of the matrix element. (see matEl_perturb_nonzero)
*/
std::complex<double> LLLhamiltonian::matEl_perturb(int i_st1,int i_st2)
{
  std::complex<double> a;
  matEl_perturb_nonzero(i_st1, i_st2, &a);
  return a;
}




/*! \fn  int LLLhamiltonian::matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a)
  
  \param i_st1, i_st2    Matrix element between the i_st1-th and i_st2-th 
                         vector of the basis. 
  \param *a              The result (matrix element).
  \returns               1 = a is nonzero, 0 = a is exactly zero (e.g. due
                         to selection rules).

  \brief This is the central function of LLLhamiltonian.
  Computes the Coulomb matrix elements between the i_st1 and i_st2-th state
   of the basis (LLL_basis). 
*/

int LLLhamiltonian::matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a)
{
  int 
    s1_diff[2], 
    s2_diff[2], 
    ind1=0, 
    ind2=0, 
    i1=0, 
    i2=0;
  
  int 
    j1 =0,
    j2 =0,
    j3 =0,
    j4 =0,
    sp1 = 0,
    sp2 = 0,
    sp3 = 0,
    sp4 = 0,
    ja = 0,
    jb = 0,
    prefactor = 0,
    parity = 0,
    nonzero=0;
  
  const int MAX_DIFF=2;   // how many 1-part. indices may be different at most

  st1=(*LLL_basis)[i_st1];
  st2=(*LLL_basis)[i_st2];
  glLogger.debug("State 1");
  glLogger.debug(*st1);
  glLogger.debug("State 2");
  glLogger.debug(*st2);
  
 
  if(!coulombME_computed) 
    {
      throw (CxErrors("No coulomb elements", __FILE__,__LINE__));
    }
  nonzero=matEl_perturb_nonzero(i_st1,i_st2,a);

  /* 
     Find out by which 1-particle states the Ne-particle states differ:
     e.g. the two 3-particle states |1 2 3> and |2 3 4> differ in one 
     state and the differing states are 1 (for the 1st) and 4 (for 
     the second). 
  */

  parity=0; // parity i1+i2+i1'+i2' (important only if st1&st2 differ at 2 pos.

  if(i_st1!=i_st2) // don't compare the states, if they are obviously identical
    {
      for(i1=0;i1<Ne;i1++){
	for(i2=0;i2<Ne;i2++){
	  //if(s1_nr(i1)+1<s2_nr(i2)) break; // a possible improvement
	  if(st1->s1_nr(i1)==st2->s2_nr(i2)) 
	    {
	      i2++;
	      break;
	    }
	}
	i2--;
	if(st1->s1_nr(i1)==st2->s2_nr(i2)) continue;
	if(ind1>MAX_DIFF-1) return nonzero;
	parity+=i1;
	s1_diff[ind1]=st1->s1_nr(i1);
	ind1++;
      }    
      for(i2=0;i2<Ne;i2++)
      {
		for(i1=0;i1<Ne;i1++)
		{
	  //if(s2_nr(i2)+1<s1_nr(i1)) break; // a possible improvement
	  		if(st1->s1_nr(i1)==st2->s2_nr(i2)) 
	  		{
	  			i1++;
	  			break;
	  		}
		}
		i1--;
		if(st1->s1_nr(i1)==st2->s2_nr(i2)) continue;
		if(ind2>MAX_DIFF-1) return nonzero;
		parity+=i2;
		s2_diff[ind2]=st2->s2_nr(i2);
		ind2++;
      } 
    }



  if(ind1!=ind2) 
    {
    	throw (CxBadValueError(__FILE__,__LINE__,"ind1 != ind2 "));
    
    }
  // Implementation of the Coulomb interaction, i.e. MAX_DIFF=2 assumed
  if(MAX_DIFF!=2) 
    {
      throw (CxBadValueError(__FILE__,__LINE__,"MAX_DIFF != 2"));
    }
  // If the two states are equal...
  if(ind1==0) 
    {
      int i,j,tmp;
      // no A_iiii term 
      // (*a)+=Ne*coulombME[0][0];
      glLogger.debug("Before adding Madelung term (%f), (%f)", 
		     a->real(), 
		     a->imag());
      (*a)+=Ne*madelung;
        glLogger.debug("Before after Madelung term (%f), (%f)", 
		       a->real(), a->imag());
      for(i=0; i<Ne; i++){
	for(j=0; j<Ne; j++)
	  {
	    tmp=abs(st1->s1_nr(i)/2-st2->s2_nr(j)/2);
	    if(i==j) 
	      {
		continue; // a^+a^+=0, aa=0
	      }
	    
	    if(((st1->s1_nr(i)-st2->s2_nr(j))%2)==0)
	      {
		/*
		  direct interaction (A_ijij) - exchange interaction (A_ijji)
		*/
		/*
		glLogger.debug(" Reading Matrixelement from Coulombmatrix at (9), (%d) = (%f), (%d), (0) = (%f)", 
			       tmp, coulombME[0][tmp], 
			       tmp, coulombME[tmp][0]);  
		*/
		(*a)+=coulombME[0][tmp]-coulombME[tmp][0];
		glLogger.debug("After adding Coulombcontributions (%f), (%f)",
			     a->real(), a->imag());	     
	      }
	    else 
	      {
	      /* 
		 if spins are different, just direct term (+assuming that
		 I will not compute matEl between |..ij> and |..ji> 
	      */
		(*a)+=coulombME[0][tmp];
	      }
	  }
      }
      //printf("->");//exit(1);
      nonzero=1;
      glLogger.debug("Diagonal Matrixelement from LLLhamiltonian (%f)",(a->real()));

      return nonzero;
    }
  // If they differ only in one one-particle st. (momentum+spin), either total 
  // momentum or total spin (z-component) cannot be equal
  if(ind1==1) 
  {
  	return nonzero;
  }
  // The states differ in two one-particle (pair) states therefore...
  j1=s1_diff[0]/2;
  j2=s1_diff[1]/2; 
  j3=s2_diff[0]/2;
  j4=s2_diff[1]/2; 
  if((j1+j2-j3-j4)% Nm)
  {
  	 return nonzero; // Total momenta must be equal
  }
  sp1=s1_diff[0]%2;sp2=s1_diff[1]%2;
  sp3=s2_diff[0]%2;sp4=s2_diff[1]%2;
  if(sp1+sp2-sp3-sp4) return nonzero;  // Total z-comp. of spin must be equal
/////////////////// Check it.
  ja=abs(j1-j4);jb=abs(j1-j3);
  prefactor=2;
  if((parity % 2)==1) prefactor=-prefactor;
  if(sp1==sp2)
    // Both pairs of one-particle states have either spins uu or spins dd
    // i.e. there's both direct and exchange term
    (*a)+=prefactor*(coulombME[jb][ja]-coulombME[ja][jb]);
  else // each of the pairs has spins ud (or du)
    if(sp1==sp3) // the two pairs are ud,ud or du,du
      (*a)+=prefactor*coulombME[jb][ja];
    else         // the two pairs are ud,du or du,ud
      (*a)+=-prefactor*coulombME[ja][jb];
// This might be improved: not all matrix elements of coulombME must be nonzero
  glLogger.debug("LLLhamiltonian: Offdiagonal element is (%f)",a->real());
  nonzero=1;                       
  return nonzero;
}


bool LLLhamiltonian::compareStates( int ist1, int ist2, int & , int &)
{
	if (ist1 == ist2)
	{
		/*
		 * Same element
		 */
	}
	else 
	{
			st1=(*LLL_basis)[ist1];
		  	st2=(*LLL_basis)[ist2];
		  	int 
		  		lIndex = 0,
		  		rIndex = 0,
		  		sign = 0;
		  		
  			glLogger.debug("State 1");
  			glLogger.debug(*st1);
  			glLogger.debug("State 2");
 			glLogger.debug(*st2);
 			st1->getSingleDiffOccNumber(*st2,rIndex,lIndex, sign);
			
	}
	return (true);
}



















