/*! \file Operator.hpp */
/*  An abstract class for representing operators and diagonalizing them.

    Implemented: diagonalization for different types of matrix
    elements (see MATRIX_TYPE).

    To be implemented: matrix element evaluation (define your own
    operator...)
*/

#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include<stdio.h>
#include<math.h>
#include<cmath>
#include<complex>
#include<limits.h>
#include<numerics/Matrix_sparse_real.hpp>
#include<numerics/matrix_full_real.hpp>
#include<numerics/matrix_full_cplx.hpp>
#include<State.hpp>
#include<Basis.hpp>
#include<myaccessories.hpp>
#include <LLLlib/Types.hpp>
#include<time.h>
#include<iostream>
#include<fstream>
 
class Operator
{

public:
  typedef enum {FULL_REAL,SPARSE_REAL,FULL_COMPLEX, SPARSE_COMPLEX} MATRIX_TYPE;

protected:
  // Expected type of the operator matrix
  

  MATRIX_TYPE matrix_type;

  
  // Stuff related to basis

  Basis *my_basis;       // Basis the op. is to be diagonalized	in
  BASIS_TYPE basisType;  // Type of basis: LC or one particular (e.g. yosBasis)
  Basis *endBasis;       // The really implemented one (yosBasis e.g.)

  FILE *logfile;         // pointer at logfile
  FILE *inpf_firstGuess; // see setFileInp_FirstGuessOfEigvecs
  FILE *outf_firstGuess; // see setFileOut_FirstGuessOfEigvecs
  

  // Stuff related to storage of matrix elements

  int matrix_allocated;          // Matrices for matEls allocated
  int spectrum_allocated;        // Matrix for spectrum allocated

  // switch to trigger saving matrices
  bool m_writeMatrixOnly;
  // Sparse matrix structures
  int *row_NrNz;         
  int *col_ind;
  double *mat_sp;
  long int nz_elements;
  double *x;
    //MaxNr of nonzero els. in sparse mat.   .  .  .   
      static const long int MAX_SIZE_SPARSE; 
  // MaxNr of Nonzero els ins parse cplx matrix
  static const long int MAX_SIZE_SPARSE_CPLX;
    // Min. value of matEl to be taken to account; set negative to 
    // use all non-zero matEls
      static const double MIN_MATEL_SPARSE; 

  // Full real matrix structures
  double *mat_full_real;
    // Max. size of the basis (refers to both FULL_REAL and FULL_COMPLEX
      static const long int MAX_DIM_FOR_FULL; 

  // Full complex matrix structures
  dComplex *mat_full_cplx;
  dComplex *x_cplx;

private:
  /*! \var static const int FREQ_REPORT_MATEL
    \brief how often should be reported
     on progress in evaluating matEl's (every x rows)
  */
  static const int FREQ_REPORT_MATEL;  

  /*! \var static const int ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE
    \brief Implicit value for ADD_TO_DIAG_WORKSPACE_SPARSE.
    \var int ADD_TO_DIAG_WORKSPACE_SPARSE
    \brief k=m+ADD_TO_DIAG_... for the SPARSE_REAL routine.
    The larger it is, the better the convergence but the more memory & 
    time per step it costs.
   */
  static const int ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE; 
  int ADD_TO_DIAG_WORKSPACE_SPARSE;

  /*! \var int ADD_TO_DIAG_WORKSPACE_SPARSE_STEP
    \brief If SPARSE_REAL fails, ADD_TO_DIAG_WORKSPACE_SPARSE will be 
    increased by this value. */
  static const int ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE_STEP; 
  int ADD_TO_DIAG_WORKSPACE_SPARSE_STEP;

  /*! \var static const char FilenameFor_tmp_eigvecsApproxs[40]
    \brief For SPARSE_REAL: if diagonalization fails (and k is to be
    increased by ADD_TO_DIAG_WORKSPACE_SPARSE_STEP), the so far
    obtained eigvec-approximations may be used starting guesses 
    for the new diagonalization. However, they first have to be
    saved into a (swap) file first and this is its filename.
    If such file cannot be opened the new diagonalization starts
    without the starting guesses.
   */
  static const char FilenameFor_tmp_eigvecsApproxs[50];


  /* If matrix elements are to be computed in a basis which comprises
     of linear combinations (see matEl_nonzero), this is used to compute
     matrix elements before they might be asked for. */
  virtual int preprocessMatEls(int newBasis_matrix_type);

  /* Private routines used by diagonalize. */
  int diagonalizeFULLREAL(int eigsToFind,double *eigvec,double *eigval); 
  int diagonalizeSPARSEREAL(int eigsToFind,double *eigvec,double *eigval); 
  int diagonalizeFULLCOMPLEX(int eigsToFind,double *eigvec,double *eigval);
  int diagonalizeSPARSECOMPLEX(int eigsToFind,double *eigvec,double *eigval); 
  
  /* See void setFileOut_FirstGuessOfEigvecs(FILE *new_outf_firstGuess) */
  void writeFirstGuessToFile(int eigsToFind,int dim,double *x);
  int getFirstGuessFromFile(int dim,double *x);
  


public:
  /* Arguments: 
     new_basis     the basis in which the op. is to be diagonalized 
     new_type      type of the operator matrix (see MATRIX_TYPE)
  */
  Operator(Basis* new_basis,int new_type);

  virtual ~Operator();

  /*!
   *
   */
  bool writeMatrixToFile(std::string& fileName, const int maxDim = -1);



  /* Set the file to which log comments should be written. Most important info is 
     also written to the screen. */
  void setLogFile(FILE *new_logfile);

  /*! \fn void setFileInp_FirstGuessOfEigvecs(FILE *inpf_firstGuess)
    \brief Diagonalization for sparse routines may be accelerated
    by giving a suitable initial guess for the eigenvectors. 
    By default this option is switched off at constructing Operator. 
    It may be switched on by calling this routine. The file is 
    expected to be of same type as files produced by 
    setFileOut_FirstGuessOfEigvecs.
   */
  void setFileInp_FirstGuessOfEigvecs(FILE *new_inpf_firstGuess);

  /*! \fn void setFileOut_FirstGuessOfEigvecs(FILE *outf_firstGuess)
    \brief If this method is called, the vectors which may be used as 
    initial guess next time will be written to file.
    See setFileInp_FirstGuessOfEigvecs.
  */
  void setFileOut_FirstGuessOfEigvecs(FILE *new_outf_firstGuess);
  
  /*! \fn change_ADD_TO_DIAG_WORKSPACE_SPARSE(int new_ADD,int new_ADD_STEP)
    \brief Change parameters ADD_TO_DIAG_WORKSPACE_SPARSE and
    ADD_TO_DIAG_WORKSPACE_SPARSE_STEP. ADD_TO... is not changed 
    if new_ADD<1, the same for new_ADD_STEP.
   */
  void change_ADD_TO_DIAG_WORKSPACE_SPARSE(int new_ADD,int new_ADD_STEP);





  /* Value of a matrix element between two states; the actual
     calculation is done by matEl_nonzero */
  std::complex<double> matEl(int,int); 

  /* Value of a matrix element between two states (returned in the 3rd
     argument); returns zero if the matEl is precisely zero. In Operator
     class it is used to display matrix elements computed in bases given 
     as linear combination of e.g. yosBasis elements. */
  virtual int matEl_nonzero(int,int,std::complex<double> *); 


  /* Guess what... Arguments:  
     eigsToFind    how many eigvals+eigvecs are to be found (note, this value may 
                   be changed according to requirements of the diagonalization 
		   routine)
     eigvec        eigenvectors (given as coeffs. of a linear
                   combination of vectors of my_basis)
     eigval        eigenvalues
  */
  virtual int diagonalize(int eigsToFind,double *eigvec,double *eigval);
};

#endif
















