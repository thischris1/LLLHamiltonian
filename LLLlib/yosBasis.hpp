/*!  \file yosBasis.hpp
  \brief See the comments in yosBasis.cpp. 
  A class for keeping a set of Yoshioka states.

   More basis-generation functions (like genBasisForConst_J_Sz) are
   anticipated for the future.
 */

#ifndef YOSBASIS_HPP
#define YOSBASIS_HPP

#include <Basis.hpp>
#include <yosState.hpp>
#include <complex>
#include <fstream>
#include <stdio.h>
#include <utils/logger.hpp>
#include <ERRORS.h>

#include <list>
#include <assert.h>
#include <iostream>
#include <vector>

class yosBasis : public Basis
{
public:
  /*! \fn  yosBasis(int new_Ne,int new_Nm,int spinYes)
   \brief   Stadard constructor (some of the basis-generating routine should
   be called after that. 
     \param new_Ne        how many electrons
     \param new_Nm        how m any flux quanta
     \param spinYes       spin polarized (0) or not necessarily polarized (1)
  */
  yosBasis(int new_Ne,int new_Nm,int spinYes);

  /*! \fn  yosBasis(FILE *in,FILE *logfile)
  \brief Get the basis from a file _with header_ (as produced e.g. by writeBasisHeaderToFile() + writeBasisToFile() ).
  */
  yosBasis(FILE *in,FILE *logfile = 0);

	yosBasis(const yosBasis &rhs);
  ~yosBasis();
yosBasis(const std::string & fileName);
  /* Construct a set of all yosStates to one total J and one total Sz
     (only for spinYes=1).
     reqJ          total J
     reqSz         total Sz
     guessDim      max. expected size of the basis;

  Output:
     determined size of the basis.
  */
  int genBasisForConst_J_Sz(int reqJ,int reqSz,long int guessDim);

 yosBasis operator = (const yosBasis &rhs);

/*! \fn int yosBasis::createAll_oneSz_oneJ(int reqJ,int reqSz,long int guessDim,yosState **output)
				  
  \brief The core of genBasisForConst_J_Sz(...) (protected method). Puts the
  generated states into output.
*/

  int createAll_oneSz_oneJ(int Ne, int Nm,
			   int reqJ,int reqSz,long int guessDim);


  int genBasisForConst_J_noSpin(int reqJ,long int guessDim);
  /* Construct a set of all yosStates to one total J
     (only for spinYes=0).
     reqJ          total J
     guessDim      max. expected size of the basis;

  Output:
     determined size of the basis.
  */

  int genBasis_allJ_noSpin();
  /* Construct a set of all yosStates to all J
     spin-polarized

  Output:
     0 -> Ok.
  */    

  int genBasis_allSz_constJ(int reqJ,long int guessDim);
  /* Construct a set of all yosStates to one total J
     with all possible spin orientations.

  Output:
     0 -> Ok.
  */    

  int genBasis_allSz_allJ();
  /* Construct a set of all yosStates with all possible spin orientations.

  Output:
     0 -> Ok.
  */    

  int setCapacity(int new_cap);
  int adjustCapacityToSize();
  int addOneState(int *new_js,int *new_spins);

	yosState * getState(long index) const;

  yosState *operator[] (long int i) const;
  //state *operator[](long int i);
  /* Returns a pointer to the i-th vector of the set. */

  /* Asking for parameters of the basis */
  int getNe() const;
  int getNm() const;
  int getSpinYes() const; 

  int increaseNm(int increment);

  /* Return the norm of the i-th state. It is automatically 1.*/
  virtual double norm(int i_st);

  /* Set this aspect ratio for all states of the basis */
  int broadcastAspectToAllStates(double new_aToB);


// Input/Output routines 
  /* File to put log information during the basis construction to */
  void setLogFile(FILE *new_logfile);  

  /*! \fn void writeBasisToFile(FILE *out)
    \brief Write down the basis into a file. Full complementary 
    to getBasisFromFile(). Some programs may however require the 
    header (use writeBasisHeaderToFile()).
  */
  void writeBasisToFile(FILE *out);

  /* Write some information about the basis into a file. */
  void writeBasisHeaderToFile(FILE *out);

  /* Add basToAdd to the basis. */
  int mergeBasis(yosBasis &basToAdd);

  /* ... in order to be able to compare it with another basis later: 
    returns 1 if the own basis and the basis in *in are the same
            0 if the bases are different (in length or in at least one state
            exits with err_msg(5) if the input file is of wrong type (not like
              a file generated by writeBasisToFile()).*/
  int amITheSameAs(FILE *in);

  /* Get the basis from a file _without header_. */
  int getBasisFromFile(FILE *in,int new_dim);

	double getAspectRatio(void) const;

  ///////////////////////////////
 // caching the wavefunctions //
///////////////////////////////

  // must be called before setWFCacheForXY
  void initWFCache();

  // initializes the cache with wavefunctions/derivatives at x,y
  void setWFCacheForXY (double x, double y);

  // set the fluxes of solenoids affecting periodic boundary-conditions
  void setWFfluxes (double alpha1, double alpha2);
 bool readFromFile(FILE * in);
  // getters for wavefunctions / derivatives
  std::complex<double> getCachedWF (int j) const;
  std::complex<double> getCachedWFdx (int j) const;
  std::complex<double> getCachedWFdy (int j) const;

  // change all j->j+shift (all Ne_ numbers in each state)
  int shiftRight_allStates(int shift);



  // Info about known quantum numbers:

  typedef struct {bool knownInAdvance; double dVal; int iVal;} Observable;

  Observable obs_totJ,obs_totSz,obs_totS;

  /*! \fn int test_totJ()
    \brief Tests whether all states of the basis have the same totJ and sets
    obs_totJ correspondingly.
   */
  int test_totJ();

  /*! \fn int test_totSz()
    \brief Tests whether all states of the basis have the same totSz and sets
    obs_totSz correspondingly.
   */
  int test_totSz();
  int createSignatures();

  inline int get_jsUpSignature(int n) {return jsUpSignature[n];}
  inline int get_jsDoSignature(int n) {return jsDoSignature[n];}
 

private:
  // Dirty thing, better don't use it (that's why it is now private)
  /* in the state i, take out the k_from spin, put it to k_to and
     shift the rest in between to the free position. */
  int shuffleSpins(int i,int k_from,int k_to);
  /* ... after shuffling, restore the original state */
  int restoreSpins(int i);
  
private:
  

  // Methods used by genBasisForConst_J_Sz(...)
  int nextSt(int *my_js);
  int nextSt(int *my_js,int Nm);

  // generate next state (without spin) recursively
  // used by genBasis_allJ_noSpin
  int nextSt_nospin(int *my_js, int tmp_Pos, int tmp_Nm);

  int totJ(int *my_js); 
  int totJ(int *my_js,int Nm); 
  //int totJ(int *my_js,int Ne, int Nm); 
  int checkSymmetry(int *spat,int *spin,int i1,int i2);
  int checkPauli_withSpin(int *spat,int *spin);

  // Methods used by genBasisForConst_J_noSpin(...)
  int checkPauli_withoutSpin(int *spat);

  // Auxiliary routine for amITheSameAs()
  int getJSCharByChar(FILE *in,int *j_l,int *spin);
  /*
   * Member functions
   */
//const long int MAX_DIM=10000;        // Maximal Nr. of basis vectors
  int Ne_;
  int Nm_;
  
  int spinYes;
  int capacity;
  
  int *jsUpSignature;
  int *jsDoSignature;        // binary shortcuts (occupation repr.) for basis states 

  FILE *logfile;

  //wavefunction-cache
  double m_WF_x;    // last x for which WF were evaluated
  double m_WF_y;    // last y ...

  std::complex<double> *m_pWaveFct;
  std::complex<double> *m_pWaveFctdx;
  std::complex<double> *m_pWaveFctdy;
  yosState **state_;
  // shuffling stuff
  int shuffled;
  yosState *st_origUnshuffled;
};



#endif














