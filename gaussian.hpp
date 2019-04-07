/*!
  \file gaussian.hpp
*/

#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP



#include <list>
#include <persist_eigSt.hpp>
#include <yosBasis.hpp>
#include <LLLhamGaussian.hpp>
#include <geometry/CxImpurity.hpp>
#include <Types.hpp>
#include <geometry/CImpurityArray.h>
const int MAX_FILENAMELENGTH = 128;

/*
 list of persist_eigSt to hold result of diagonalisation
*/
typedef std::list<persist_eigSt> TRESULT;

/*!
  \class gaussian
  \brief Class controlling calculation of eigenstates in a fqh-system containing
  a gaussian. Reading of parameter file 'gaussian.par', writing of files with
  results of calculation: basis, eigenvectors and -values (names according to gaussian.par).
*/
class gaussian {

private:
  /*!
    \const MAX_FILENAMELENGTH
    \brief Maximum length of filenames expected to be in 'gaussian.par'.
    Constraint: >= 30 because of appended parameters !
  */
	bool m_potentialOnly;
  //const static int MAX_FILENAMELENGTH;

public:
  /*!
    \fn gaussian ()
    \brief Standard constructor
  */
  gaussian (bool n_potentialOnly);

  /*!
    �\fn int readParameters (char *name_of_paramfile)
    \brief Reads parameter file.
    \param name_of_paramfile Name of the ascii file containing information about the system 
    (referred to as 'gaussian.par').    
    \return 0 if successful


    Format of Parameter file is as follows:
    "../out/bs"  # root-name of basis file
    "../out/vc"  # root-name of vector file
    "../out/dn"  # root-name of density file
    "../out/ld"  # root-name of landau-diagonal file
    "../out/pot"	# root-name of potential file
    3		# Ne: Nr. of electrons
    9		# Nm: Nr. of flux quanta (i.e. Ne/Nm=filling factor)
    0 		# spinYes: 0=spin polarized, 1=not necessarily sp. pol.
    0		# reqSz: dtto, with total Sz (applies only if spinYes=1)
    0		# mat_type: 0=FULL_REAL, 1=SPARSE_REAL, 2=FULL_CPLX
    1		# type of vector-file to generate: 0->ascii, 1->raw binary, 2->fortran's unformatted
    10		# eigsToFind: Nr. of eigvals/eigvecs to be found
    1.0		# a: size of the system (vert.)
    1.0		# b: size of the system (horiz.)
    0.0		# bli: related to finite thickness
    0		# type of gaussian potential: 0 -> gaussian, 1 -> delta
    0		# type of e-e interaction: 0 -> Coulomb, 1 -> hardcore
    0.5		# Strength of gaussian
    0.1		# Width of gaussian
    0.5		# x-Position of gaussian
    -0.5		# Strength of hole
    0.2		# Width of hole
    0.5		# y-Position of hole
    0.0		# energy-offset
    100		# xkmax: Sum from -kmax to kmax for Gaussian in x-direction (resp. hole)
    100		# ykmax: Sum from -kmax to kmax for Gaussian in in y-direction (resp. hole)
  */
  int readParameters (std::string & name_of_paramfile);

  /*!
    \fn void createStandardFilenames (char *strExtraAppendix)
    \brief Creates Filenames out of read in root-names and parameters
    \param strExtraAppendix will be appended directly before file suffix.
  */
  void createStandardFilenames (std::string & strExtraAppendix);

  /*!
    \fn int calculate ()
    \brief Perform the calculation for the current parameters of the system.
    Saves results in files directed by the root names in 'gaussian.par' using format also specified there.
  */
  int calculate ();

  /*!
    \fn int writePotential ()
    \brief Writes a file containing the shape of the used gaussian-potential in ascii-format
    like x, y, V(x,y) (gnuplot readable)
    The filename will be constructed according to the root name given in 'gaussian.par'
    \return 0 if successful
  */
  int writePotential ();

  /*!
    \fn int writeEvalDensity ()
    \brief Creates a file 'evalDensity.par' that contains the parameters to evaluate
    the density and cutrrent-density of the groundstate.
    \return 0 if successful
  */
  int writeEvalDensity ();

private:
  /*
    Private members
  */

  /*!
    \fn int saveVectors (TRESULT EigVectors)     
    \brief saves eigenvalues and -vetors in file.
  */
  int saveVectors (TRESULT & EigVectors);

  /*!
    \fn int diagonalize (TRESULT listEigVectors, yosBasis & basis)
    \brief Will be called by calculate to diagonalize the Hamiltonian.
    \param listEigVectors list of persist_eigSt to store found eigenvectors in
    \param basis Basis to calculate in
    \return number of found eigenvectors
  */
  int diagonalize (TRESULT & listEigVectors, yosBasis & basis);
    
private:  
  /*
    Data
  */


  /*
    Parameters controlling the basis: 
    m_iNe = number of electrons
    m_iNm = number of flux quanta
    m_iSpinYes = Spin 1->yes, 0->no
    m_iReqSz = which z-component of spin is required (meaningful only if spinYes = 1)
    m_bComplex = true->use complex eigenvectors
   */
  int 
    m_iNe,
    m_iNm,
    m_iSpinYes,
    m_iReqSz;

  bool
    m_bComplex,
    isRandomArray,
    isRandomMatrix; //! Use a Random Landau Matrix (RLM)
/*
 * Possible parameters for  RLM
 */
 	double m_impStrength;
 	double m_correlationLength;
	int m_randomIMax;
  /*
    Paramters for LLLhamGaussian:
    m_iMatType = type of matrix to be used
    m_iPotentialType = 0->gaussian, 1->delta - gaussian
    m_iInteractionType = 0->coulomb, 1->hardcore
    m_iXKMax = maximum k to sum over in direction of x for calculating matrix-elements
    m_iYKMax = maximum k to sum over in direction of y for calculating matrix-elements
  */  
  int   
    m_iMatType, 
    m_iPotentialType,
    m_iInteractionType,
    m_iXKmax,
    m_iYKmax;


  double 
    m_dA, 
    m_dB, 
    m_dBli,
    m_dEnergyOffset,
    m_dAlpha1,
    m_dAlpha2;
    
    CImpurityArray m_impArray;


  /*
    other parameters
    m_iVectorFileFormat of Vector's file. 0->ascii, 1->binary, 2->fortran's unformatted
    m_iEigsToFind number of eigenvalues/vectors to find
  */
  int 
    m_iVectorFileFormat,
    m_iEigsToFind;
 
  /*
    Filenames in order of apperance:
    Basis
    Eigenvectors
    Density (used by writeEvalDensity)
    Potential (used by writePotential)
    Occupation-numbers (Landau-matrix, used by writeEvalDensity)
  */
  std::string  m_strFileBas;
  std::string  m_strFileVec;
  std::string  m_strFileDens;
  std::string  m_strFileV;
  std::string  m_strFileLand;
  std::string  stateFilePrefix;
  std::string m_spectrumFilePrefix;
  std::string  baseFileprefix;
  std::string  densFilePrefix;
  std::string  landauFilePrefix;
};

#endif
