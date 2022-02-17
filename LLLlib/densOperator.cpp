/*!
 \file densOperator.cpp
 \brief Implementation of density operator based on OneParticleOperator. Provides particle density at
 arbitrary points inside unit cell and occupation numbers of one-particle states.
*/

#include <densOperator.hpp>
#include <utils/logger.hpp>

DensOperator::DensOperator ( const Basis &basis_to_use, 
		 const eigSt &state_to_eval, 
		 int mat_type ) :
  OneParticleOperator (basis_to_use, state_to_eval, mat_type, true, false)
{
}
 
DensOperator::DensOperator ( const Basis &basis_to_use, 
		 const eigSt &state_to_eval, 
		 int mat_type,
		 bool bDoNotComputeLandauMatrix) :
  OneParticleOperator (basis_to_use, state_to_eval, mat_type, true, bDoNotComputeLandauMatrix)
{
}
 
DensOperator::DensOperator ( const Basis &basis_to_use, 
		const OneParticleOperator &ancestor, 
		int mat_type ) :
  OneParticleOperator (basis_to_use, ancestor, mat_type, true)
{
}


DensOperator::~DensOperator ()
{
}

std::complex<double> DensOperator::getOneParticleEl (int i, int j)
{
  return conj( endYosBasis->getCachedWF(i) ) * endYosBasis->getCachedWF(j);
}


int DensOperator::getOccupationNumbers (TOCCNUMLIST & occnumlist)
{

  int ifail = -1;
  int lcwork = 16*Nm;  // 64 would give optimal performance on most computers
  int lrwork = 3*Nm;  // at least 3*Nm
  

  std::complex<double> *tmpDens = new std::complex<double>[Nm*Nm];
  double *tmpEigs = new double[Nm];
  double *rwork = new double [lrwork];
  std::complex<double> *cwork = new std::complex<double>[lcwork];

  assert (tmpDens != NULL);
  assert (tmpEigs != NULL);
  assert (cwork != NULL);
  assert (rwork != NULL);

  glLogger.info ("Diagonalizing density-matrix...");


  // implicitly assumed...
  assert (2*sizeof (double) == sizeof (std::complex<double>));



  memcpy (tmpDens, getLandauMatrix(), Nm*Nm*2*sizeof(double));
  /*
  for (int i=0; i<Nm; i++)
    {
      for (int j=0; j<=i; j++)
      {
	std::cout << tmpDens [i*Nm+j] << ' ';
      }
      std::cout << '\n';
    }
  */


  // real, symmetric
  //f02faf_ (&job, &uplo, &Nm, tmpDens, &Nm, tmpEigs, work, &lwork, &ifail);

  // complex, hermitian

  glLogger.info ("done, ifail = %d", ifail);
  
  if (ifail == 0)  // successsul
    {
      TOCCNUM occnum;  
      double maxcoef;
      int maxj=-1;
      double tmp;
      
      for ( int i=0; i<Nm; i++ ) 
	{
	  
	  // find maximal coefficient in eigenvector
	  maxcoef = 0.0;
	  for ( int k=0; k<Nm; k++ )
	    if ( ( tmp = abs(tmpDens[i*Nm + k]) ) > maxcoef )
	      {
		maxcoef = tmp;
		maxj = k;
	      }
	  
	  // storing results in the list
	  occnum.maxj = maxj;
	  occnum.occ = tmpEigs[i];
	  occnum.diag = getLandauMatrix()[(Nm+1)*i].real();
	  occnumlist.push_back( occnum );
	
      }
      
      // sort occupation numbers by OP-WF to which they contribute most
      //occnumlist.sort();
      
    } // of if (ifail...)
 
  delete[] tmpDens;
  delete[] tmpEigs;
  delete[] rwork;
  delete[] cwork;

	
  return ifail;
}

bool operator< (const TOCCNUM & occ1, const TOCCNUM & occ2) 
{ 
  return (occ1.maxj < occ2.maxj); 
}



















