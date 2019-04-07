/*! \file LLLhamiltonian.hpp */
/* A class for representing Hamiltonians of the lowest Landau
   level. my_basis (see Operator.hpp) must be a many-spin-particle Yoshioka
   basis (rectangular geometry, periodic BC).
   
   Implemented: evaluation of matrix els of the Coulomb interaction and hardcore interaction
   
   Not implemented: matrix els of an inhomogeneity (perturbation) if needed   
 */

#ifndef LLLhamiltonian_HPP
#define LLLhamiltonian_HPP

#include<Basis.hpp>
#include<yosBasis.hpp>
#include<Operator.hpp>
#ifdef NAG
// NAG fortran library routine for incomplete Gamma function
extern "C" void s14baf_(double *,double *,double *,double *,double *,int *);
#endif

/*! \class LLLhamiltonian : public Operator
  \brief A class for representing Hamiltonians of the lowest Landau
   level. my_basis (see Operator.hpp) must be a many-spin-particle Yoshioka
   basis (rectangular geometry, periodic BC).
   
   Implemented: evaluation of matrix els of the Coulomb interaction and hardcore interaction
   
   Not implemented: matrix els of an inhomogeneity (perturbation) if needed   
 */
class LLLhamiltonian : public Operator
{
protected:

  /* Physics-related:
   */
  const static double CUTOFF_MIN_Q2;   // in computeCoulomb: cutoffs for
  const static double CUTOFF_MAX_Q2;   // q=(q_x,q_y) (in MIN<|q|^2<MAX)

  yosBasis *LLL_basis;                // basis of the space
  int Ne, Nm, spinPolarized;          /* Nr. of electrons, Nr. of flux quanta
					 (determined from LLL_basis) */
  double bli;                         // ~thickness of the 2DEG
  double a,b;                         // size of the system (rectangle)
 
  double madelung;                    // diagonal term (W, Chakr. p44)

  /* Numerics-related:
   */
  int coulombME_computed;           // have the matrix coulombME been compd?
  double **coulombME;
  double computeMadelung();

  /* Related to evaluation of matrix elements 
   */
  const yosState *st1,*st2;          
  /* 1-electron states are ordered by 2*j+spin (spin=0,1). */
  /*
  int s1_nr(int i) const; 
  int s2_nr(int i) const;
*/

public:
  /*!
    \fn LLLhamiltonian(yosBasis *new_basis,int new_type,double bli_new,
		 double a_new, double b_new, bool bHardcore = false);
		 
    \param new_basis                    basis in which the matEls are to be computed
    \param Ne_new, Nm_new, spin_yes     description of the Yoshioka basis; it is maybe
                                      a good idea to replace it by a pointer to basis
    \param bli_new,a_new,b_new          (bli) ~thickness of the 2DEG, (a,b) size of the rectangle
    \param bHardcore                   wether or not to use hardcore-interaction (false->coulomb)
  */
  LLLhamiltonian(yosBasis *new_basis,int new_type,double bli_new,
		 double a_new, double b_new);

  virtual ~LLLhamiltonian();

  /* Comp. the matEls of the 2-particle Coulomb interaction operator */
  int computeCoulomb();

  /*!
    \fn int computeHardcore()
    \brief Computes two-particle matrix elements for hardcore interaction.
  */
  int computeHardcore();

  /*!
   * \fn computeNoInteraction()
   * \brief just fill the appropriate array and set computeColoumb to 0
   */
  int computeNoInteraction();
  /* Coulomb matrix elements (overloaded Operator::matEl_nonzero(...)) 
   */
  int matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a);

  /* Matrix element of inhomogeneity (perturbation); set to zero here,
     may be redefined (see also Operator::matEl(State *st1,State *st2)).
   */
  virtual int matEl_perturb_nonzero(int,int,std::complex<double> *);
  std::complex<double> matEl_perturb(int,int);
  //! compare 2 states from the basis 
  bool compareStates (int ist1, int ist2, int & ind_state1, int & ind_state2);
  double get_Lx(void) const {return (a);};
  double get_Ly(void) const {return (b);};
};




/* Related to evaluation of matrix elements 
 */

#endif











