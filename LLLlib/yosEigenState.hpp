#ifndef YOSEIGENSTATE_HPP
#define YOSEIGENSTATE_HPP

/*!
	\file        yosEigenState.hpp
	\brief The file contains the declaration of the class yosEigenState
	\author     Christian Mueller
	\date       07 May 03

	
*/ 
#include <LLLlib/Types.hpp>
#include <mkl_rci.h>

 
#include <vector>
#include <LLLlib/persist_eigSt.hpp>
#include <LLLlib/TPositionArray.h>
#include <LLLlib/LLLlib.h>

#include <LLLlib/yosBasis.hpp>
#include <LLLlib/Types.hpp>
#include <geometry/CxPosition.hpp>
#include <geometry/CxPeriodicPosition.hpp>
#include <complex>
#include <string>
#include <fstream>
#include <geometry/CxVortex.hpp>
#include <numerics/matrix_full_cplx.hpp>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>




extern "C" 
{
  dComplex reducedWfct(int, int, double, double);
  /*!
    NAG function to search minima
  */

#ifdef NAG
void e04usf_(int *, int *, int *, int *, int *, int*, int *, int *, double *, double *,  double *, double *, 
		void (*pCon)(void),
		void (*pObj)(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *),
		int *, int *, double *, double *, double *, double *, double *, double *, double *,
	        double *, int *, int *, double *, int *,int*, double *, int *);
 


  void e04unf_(int *, int *, int *, int *, int *, int*, int *, int *, double *, double *,  double *, double *, 
		void (*pCon)(void),
		void (*pObj)(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *),
		int *, int *, double *, double *, double *, double *, double *, double *, double *,
	        double *, int *, int *, double *, int *,int*, double *, int *);
  /*!
    NAG dummy function
  */
  typedef void (*DUMMYPTR)();
  /*!
    Methods needed for communication with the NAG library
   */
  void e04udm_();

  void e04urf_(char *, int);
// specify output here
   void x04aaf_(int*, int *);
  
#endif //NAG
  
  typedef void (*OBJPTR) (int *, int *, int *,  int *, double *, double *, double *, int *, int *, double *);
  typedef void (*CONFPTR)(int *, int *, int *, int *, int *, double *, double *, double *);
  #ifdef GSL
    double my_f (const gsl_vector *v, void *params);
    void my_df_gsl (const gsl_vector * x, void * params, gsl_vector * g);
    void my_df(double *x, void *params, double * fjac);
    void my_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g);
  //	void * parameter = 0;
#endif //GSL
  #ifdef MKL
		double my_f (CxPosition &);
		void my_df(double *x, void *params, double * fjac);
		void my_df_gsl (const gsl_vector * x, void * params, gsl_vector * g);
#endif    // MKL 

}
void setDoubleArray (double *, int, double);



/*!
  \brief convenience function for comparing complex numbers
*/
bool complexSmaller( std::complex<double> , std::complex<double>);
/*!
  \brief  convenience function for comparing complex numbers
*/
bool complexLarger( std::complex<double> , std::complex<double>);



/*!
  \class yosEigenState : public persist_eigSt
  \brief An eigenstate with a yosBasis
  Allows to find vertices and many more (to come)
*/

class yosEigenState : public persist_eigSt 
{
public:
  //! Constructor 
  yosEigenState(yosBasis *, std::complex<double> *);
  
  //! Constructor which uses the output of Operator::diagonalize
  yosEigenState(yosBasis *newBasis, double * eigenvec, 
		double * eigenVal, int stateNo);
  
  //! Constructor from persist_eigSt
  yosEigenState ( yosBasis *, persist_eigSt &);
  //! Constructor with explicitly set fileName
  yosEigenState ( yosBasis *, persist_eigSt &, std::string);
  //! Ctor whith electron positions explicitly passed
  yosEigenState (yosBasis *, persist_eigSt &, std::vector<CxPeriodicPosition> &);
	//! Constructs with a file name for the state (open from file)
  yosEigenState (yosBasis *, std::string stateFile, int stateNo=0);
	//! copy ctor  
  yosEigenState( const yosEigenState &rhs);
  //! gets the next state from stream
 yosEigenState( yosBasis *, std::istream & inStream);
  //! Destructor
  virtual ~yosEigenState();
  //! Adds Electron position (Ne-1) to state
  void addElectronPositions(CxPeriodicPosition * );
  //! Reads electronPosition from File
  void readElectronPositions (std::string);
  //!Deletes all ElectronPositions currently available
  void resetElectronPositions();
  //! Resets all immediate values.
  bool reset(void);
  //! Resets all electronpositions
  bool setElectronPositions (int, CxPeriodicPosition *);
  //! Sets sample electrons to values passed in vector
  bool setElectronPositions (std::vector<CxPeriodicPosition> &);
  //! Randomly generates electron positions
  bool generateRandomPositions(double & );
  //! Shall be able to calculate the windingnumber at a given point
  double getWindingNumber(CxPeriodicPosition &) const;
  
  //! Move the electrons in a star like way (multiply the vectors)
  void moveElectrons(const double factor, const CxPeriodicPosition &);
  //! Does the actual work for us
  void findZeros(void);
 
  //! Gets the position of the zeros found
  std::vector<CxPeriodicPosition> getShallowZeros(void) const;
  //! Gets vortices
  std::vector<CxVortex> getShallowVortices(void)const;
 //! Writes  the quasi-sp wavefunction into a file
  int writeReducedWfct(int, std::string , double ) const;
  //! Methods to calculate the reduced sp-wavefunctions
  std::complex<double> sumofSPWavefunctions(const int , CxPosition & ) const;
  //! Convenience Method 
  std::complex<double> sumofSPWavefunctions(const int , double , double ) const;

  std::complex<double> sumofSPWavefunctions(const int , double[] ) const;
  std::complex<double> sumofdySPWavefunctions(const int , CxPosition & ) const;

  std::complex<double> sumofdxSPWavefunctions(const int , CxPosition & ) const;
  //! Calculates the sub-determinants for the quasi sp-wavefunction
  bool calculateFixedPart(void);
  //! Gets the aspect ratio of the underlying system
  double getAspectRatio(void) const;
  //! returns the zeros as a stdlib vector
  std::vector <CxPeriodicPosition> getVectorOfZeros(void) const;
  //!
  std::vector <CxPeriodicPosition> getElectronPositions(void) const;
  //! Gets the number of electrons in the state (from the basis)
  int getNe(void)const;

  //! returns the weight of the wavefunctions when fixing electron positions.
  float getWeight(void);
	void initializeSearch(void);
static yosEigenState * getInstance(); 
	bool readState(std::istream & in);
protected:
  

  double sampleElectronSite(CxPeriodicPosition & ) const;
  //! Later possibly returning the zero-positions
  
  int searchAllZeros (void);

  //! Use the NAG librray to do the searching
#ifdef NAG
	void initializeNag();
  	int searchNag(void);
#endif
#ifdef INTEL
		
#endif
	void initializeIntel(void);
		int searchIntel(void);
	bool freezeInstance ();
	
	int searchGsl(void);
private:
  
  yosEigenState();
  //! Helper initializer
  bool init (yosBasis *, std::complex<double> *);
  //! Basis 
  yosBasis *myBasis;
  //! Array of electronSites 

  //! electronSites in  STL vector
  std::vector<CxPeriodicPosition> electronVector;
  //! Array of zeros found

  TPositionArray <CxPeriodicPosition> zeroPositions;

    TPositionArray <CxVortex> vortexPositions;

  int noOfElectrons;
//  int numberOfZeros;
//  int numberOfVortices;
  //!
  double aspectRatio;
  //!
  std::complex<double> *fixedPart; 
  bool initialized;
  static const int curvIntStepsMax;
  static const double curvRadius;
  static const double vortexTolerance;
  //! this is a trial now
  static yosEigenState *m_instance;
   bool isSearchInitialized;
   #ifdef MKL
   // Searcher variables
   _TRNSPBC_HANDLE_t m_solverHandle;
	
   #endif
   double *fvec; //! for intel solver
	double *fjac; //! jacobi matrix for intel solver
	double eps[6];  //! for jacobi matriux
};

  //////////////////////
 // friend operators //
//////////////////////

  // stream inserter
  std::ostream & operator<< (std::ostream &out, const yosEigenState &eigst);

  // stream extracter
  std::istream & operator>> (std::istream &in, yosEigenState &eigst);



#endif 
