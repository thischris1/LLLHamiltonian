/*!
 \file densOperator.hpp
 \brief Declaration of density operator based on OneParticleOperator. Provides particle density at
 arbitrary points inside unit cell and occupation numbers of one-particle states.
*/


#ifndef DENSOP_HPP
#define DENSOP_HPP

#include <OneParticleOperator.hpp>
#include <cassert>
#include <list>


/*!
  \typedef TOCCNUM
  \brief structure holding occupation-numbers occ after diagonalization of density-matrix
  and the quantum number maxj of the one-particle-wavefunction this number contributes most to.
  diag is the diagonal element of the landau matrix before diagonalizing.
*/
typedef struct
{
  /*!
    quantum number j
  */
  int maxj;
  /*!
    eigenvalue, to which maxj contributes most
  */
  double occ;
  /*!
    diagonal element dens(maxj, maxj) of density-matrix
  */
  double diag;
}TOCCNUM;

/*!
  
  compares two structs of type TOCCNUM by maxj; needed to sort TOCCNUMLIST
*/
bool operator< (const TOCCNUM & occ1, const TOCCNUM & occ2);

/*!
  \typedef TOCCNUMLIST
  \brief list of occupation numbers after diagonalization of density-matrix
*/
typedef std::list<TOCCNUM> TOCCNUMLIST;

#ifdef NAG
// diagonalizing real symmetric matrix
extern "C" void f02faf_ (char *, 
			 char *, 
			 int *, 
			 double *, 
			 int *, 
			 double *, 
			 double *, 
			 int *, 
			 int *);

// diagonalizing complex hermitian matrix
extern "C" void f02haf_ (char *, 
			 char *, 
			 int *, 
			 std::complex<double> *, 
			 int *, 
			 double *, 
			 double *,
			 std::complex<double> *, 
			 int *, 
			 int *);

#endif
/*!
 \class DensOperator : public OneParticleOperator
 \brief class representing density-operator that can be evaluated at distinct points
 inside the unit cell and which can be diagonalized to obtain occupation numbers \see getOccupationNumbers
*/
class DensOperator : public OneParticleOperator
{

private:


public:

  /*! \fn DensOperator ( const Basis &basis_to_use, 
	   const eigSt &state_to_eval, 
	   int mat_type )
      \brief
    basis_to_use: Basis to be used
    state_to_eval: State, for which density is to be evaluated
    mat_type: type of density matrix (only 0 = full_real implemented)
  */
  DensOperator ( const Basis &basis_to_use, 
	   const eigSt &state_to_eval, 
	   int mat_type );


  /*! \fn DensOperator ( const Basis &basis_to_use, 
	   const eigSt &state_to_eval, 
	   int mat_type,bool bDoNotComputeLandauMatrix )
      \brief
    The same as for the ctor above, intended for use in derived classes,
    which should compute their own Landau matrix in their ctor.
  */
  DensOperator ( const Basis &basis_to_use, 
	   const eigSt &state_to_eval, 
	   int mat_type,bool bDoNotComputeLandauMatrix );

  /*! \fn DensOperator ( const Basis &basis_to_use, 
	   const OneParticleOperator &ancestor, 
	   int mat_type)
    \brief
    basis_to_use: Basis to be used
    ancector: another OP-operator from which the landaumatrix can be inherited
    mat_type: type of density matrix (only 0 = full_real implemented)
  */
  DensOperator ( const Basis &basis_to_use, 
	   const OneParticleOperator &ancestor, 
	   int mat_type);

  /*!
    Destructor
  */
  ~DensOperator ();

  /*!
    One-particle matrix element of this operator
  */
  std::complex<double> getOneParticleEl (int i, int j);

  /*!
    Calculates occupation numbers by diagonalizing the landau-matrix
    returns 0 if successful
  */
  int getOccupationNumbers (TOCCNUMLIST & occnum);

};

#endif















