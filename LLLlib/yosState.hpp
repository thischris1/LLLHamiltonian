/*!
  \file yosState.hpp
  \brief This file holds the declaration of class yosState
  \author Karel Vyborny
*/

#ifndef YOSSTATE_HPP
#define YOSSTATE_HPP

#include<State.hpp>
#include <stdio.h>
#include <complex>


//! k Sum in the SP-Wavefunctions k = (-kSPsum...kSPsum)
const int kSPsum = 10;


/*! 
  \class  yosState : public State
  \brief A class for keeping one few-electron Yoshioka state with or without spin. 

*/
struct cmp_yosStates_differences{
  static const int MAX_DIFF = 2;
  int j_diff_l   [MAX_DIFF];
  int spin_diff_l[MAX_DIFF];
  int pos_diff_l [MAX_DIFF];
  int j_diff_r   [MAX_DIFF];
  int spin_diff_r[MAX_DIFF];
  int pos_diff_r [MAX_DIFF];};
class yosState : public State
{
public:

  /* Initialize a state without spin.
     new_Ne    Nr. of electrons
     new_Nm    Nr. of flux quanta
     new_j     j's of the electrons
   */
  yosState(int new_Ne,int new_Nm,int *new_j);

  /* Initialize a state with spin.
     new_Ne    Nr. of electrons
     new_Nm    Nr. of flux quanta
     new_j     j's of the electrons
     new_spin  spins of the electrons
   */
  yosState(int new_Ne,int new_Nm,int *new_j,int *new_spin);

  virtual ~yosState();
  //! Copy constructor
  yosState(const yosState &);

  //! Copying via copy ctor
  yosState & operator = (const yosState &rhs);



  /* Returns j of the i-th electron. */
  int j(int i) const {return j_[indInRange(i)];}

  //! Returns spin of the i-th electron. By convention 1=up, 0=down.
  int spin(int i) const {if(spinYes) return spin_[indInRange(i)]; else return 0;}

  //! Is the state spin polarized? (yes = it has no spin) 
  int hasSpin() const {return (int) spinYes;}

  //! Returns Number of electrons Ne. 
  int getNe() const {return (int) Ne_;}

  //! Returns Number of flux quanta Nm. 
  int getNm() const {return (int) Nm_;}
//! Increases Nm_.
  int increaseNm(int increment);
  //! Sets the aspect ratio (influences only SPwaveFct); 
  int setAspect(double new_aToB);

  //! Returns the aspect ratio of the system/wavefunction

  double getAspect(void) const { return aToB;};

  /* Print the state in a symbolic form |......> */
   void print() const;
 // Text output; basic routines: printing to a string
  void sPrint_num(char str[],int len) const; 
  void sPrint_occ(char str[],int len) const;

  /* Print the state in a symbolic form |......> to a file */
  void fprint(FILE *f) const;

 /* Print the state in a symbolic form |1+2-2+5+> */
 
  /* Print the state in a symbolic form |u d b> to a file */
  virtual  void fprint_x(FILE *f) const;

  /*!
    \fn void setSolenoidFluxes (double alpha1, double alpha2)
    \param alpha1 flux of solenoid parrallel to y-axes in units of hbar/e
    \param alpha2 flux of solenoid parrallel to x-axes in units of hbar/e
    \brief Changes the fluxes of solonoids in the torus to modify the periodic
    boundary conditions. (see paper of Tao & Haldane: topological invariant & impurity-effect)
    This affects the one-particle wavefunctions and their derivatives.
    Default values are set in constructor to 0,0.
  */
  void setSolenoidFluxes (double alpha1, double alpha2);

  double getSolenoidFlux1() const;
  double getSolenoidFlux2() const;
  

  /* One particle wavefunction: value of the j-th function at the place (x,y)
   */
  std::complex<double> SPwaveFct(int j,double x,double y) const;

  /* One particle wavefunction: derivative with respect to x of the j-th function at the place (x,y)
   */
  std::complex<double> SPwaveFctdx(int j,double x,double y) const;

  /* One particle wavefunction: derivative with respect to y of the j-th function at the place (x,y)
   */
  std::complex<double> SPwaveFctdy(int j,double x,double y) const;

  /*!
    \fn compare(const yosState) const;
    compare this state with another one
    \param rhs the other one
    \return true if equal
  */
  virtual bool compare (const yosState& rhs) const;
  
  
  virtual int compareAndFindDiffs(const yosState &rhs,
				  cmp_yosStates_differences &diff) const;
  
  /*! Makes a deep copy of the occ states, thus detaching the object. After detaching you can independently handle this object from the older ones. Example
  yosState a(1,2,3);
  yosStae b;
  b = a;
  // If a dies (runs out of scope b does not have any access to the occ numbers. To prevent this:
  b.detach();
  \return false if any error occurs
  */
  bool detach(void);
//! needed for matrix evaluation
int s1_nr(int i) const {return (int) j(i)*2+spin(i);};
// needed for matrix evaluation
int s2_nr(int i) const {  return (int) j(i)*2+spin(i);};

  /*!

  \fn bool operator == (const yosState&) const
  \brief overloaded comparison operator 
  \todo put this into State?
  \param rhs the other one
  \return true if equal
  */
  
  bool operator == (const yosState& rhs) const;

  /*!

  \fn int getDifferentNumbers( const yosState &rhs, int& rIndex, int& lIndex, int& sign  ) const
  
  */
  int getDifferentNumbers(const yosState &rhs, 
			  int& rIndex, 
			  int& lIndex,
			  int& sign  ) const;  
/*!
    
  \fn bool getSingleDiffOccNumber ( const yosState &rhs,  int& rIndex,  int& lIndex,  int& sign  ) const;
  */
  bool getSingleDiffOccNumber ( const yosState &rhs,  
				 int& rIndex,  
				 int& lIndex,  
				 int& sign  ) const;
  

  /* For safe access to the j[] and spin[] arrays. Returns i if 
     i\in\{0,..,Ne_-1\} and zero otherwise. */
  int indInRange(int i) const;

/*! \fn  void shiftRight_State(int shift)

\brief Shift all J's of the state by shift to right (mod Nm_);
 */
  void shiftRight_State(int shift);

/*! \fn int st_nr(int i)
  \brief Attributes a unique number to the i-th one particle state.
*/
  int st_nr(int i) const;

protected:
  int * getjStates(void) const ;
  int * getSpinStates(void) const;


private:
  int Ne_;
  int Nm_;
  int spinYes;
  int *j_;
  int *spin_;
  double aToB;         // aspect ratio
 // static int reference; //! Needed for copying
  double alpha1_;     // flux in units of hbar/e of solenoid parallel to y-axis
  double alpha2_;     // flux in units of hbar/e of solenoid parallel to x-axis
};







  /*!    
    \fn  int equalSpat(const yosState &st1,const yosState &st2);
  */
  /* Compares the spatial parts of the two yosStates. Returns
     1 if they are equal and returns 0 if not (or if they are even
     not of the same type). */
int equalSpat(const yosState *st1,const yosState *st2);


  /*!    
    \fn not-a-friend int equalType(const yosState &st1,const yosState &st2);
  */
  /* Compares two yosState's st1, st2. Returns 1 if Ne's are
     the same and Nm's are the same. Returns 0 otherwise.*/
int equalType(const yosState *st1,const yosState *st2);






#endif












