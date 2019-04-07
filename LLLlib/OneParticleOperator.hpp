/*!
  \file OneParticleOperator.hpp
*/
#ifndef ONEPARTICLEOPERATOR_HPP
#define ONEPARTICLEOPERATOR_HPP

#include <Operator.hpp>
#include <Basis.hpp>
#include <yosBasis.hpp>
#include <eigSt.hpp>
#include <cassert>

/*!
 \class OneParticleOperator : public Operator
 \brief class to calculate expectation values of one-patricle-operators in a multi-particle state
 encapsulates clacuation of 'landau-matrix'
*/
class OneParticleOperator : public Operator
{

public:
  
  /*!
    \fn OneParticleOperator(const Basis &newBasis, const eigSt &state, int new_type, bool bNeedsWF = false, bool bDoNotComputeLandauMatrix = false)
    \brief newBasis: basis to construct Operator for
    state: state in which operator will be evaluated
    new_type: how to save matrix (full, sparse)
    bNeedsWF: whether the oparator needs the wavefunctions/derivatives to be computed
  */
  OneParticleOperator(const Basis &newBasis, const eigSt &state, int new_type, bool bNeedsWF = false, bool bDoNotComputeLandauMatrix = false);

  /*!
    \fn OneParticleOperator(const Basis &newBasis, const OneParticleOperator &ancestor, int new_type, bool bNeedsWF = false);
    \brief newBasis: basis to construct Operator for
    ancestor: OP-operator constructed for the same state to inherit landau-matrix
    new_type: how to save matrix (full, sparse)
    bNeedsWF: whether the oparator needs the wavefunctions/derivatives to be computed
  */
  OneParticleOperator(const Basis &newBasis, const OneParticleOperator &ancestor, int new_type, bool bNeedsWF = false);

  virtual ~OneParticleOperator();
  
  /*!
    \fn virtual int setXY(double new_x,double new_y)
    \brief set position for evaluation
  */
  virtual int setXY(double new_x,double new_y);
  
  /*!
    \fn double evalAt (double x, double y)
    \brief evaluates Operator at Position x,y
  */
  double evalAt (double x, double y);

  /*!
    \fn virtual int matEl_nonzero(int i_st1,int i_st2, std::complex<double> *a)
    \brief must be overridden by derived class to implement specific operator
  */
  virtual int matEl_nonzero(int i_st1,int i_st2, std::complex<double> *a);
  
  virtual int preprocessMatEls(int newBasis_matrix_type);

  /*!
    \fn const double *getLandauMatrix () const
    \brief retrieves a constant pointer to the Landau-Matrix
    is needed to share landaumatrices when object is created by second constructor 
  */
  const std::complex<double> *getLandauMatrix () const { return m_pCmplxLandauMatrix; };

protected:

  /*!
    \fn virtual std::complex<double> getOneParticleEl (int i, int j)
    \brief must be overridden by derived class to implement specific operator
    i,j quantum number of one-particle states
    \return allways 0
  */
  virtual std::complex<double> getOneParticleEl (int i, int j) { return  std::complex<double>(0.0,0.0); };

  // constructs ALL-Landaumatrix
  int setupLandauMatrix(void);
  
  /*!
    \fn int setupLandauMatrixForState (void)
    \brief constructs Landau-Matrix for the state, the operator was created for
  */
  int setupLandauMatrixForState (void);

  /*!
    \fn virtual double getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2)
    \brief retrieved matrixelement of operator
    k1: ?
    k2: ?
    i1: ?
    i2: ?
  */
  virtual std::complex<double> getOneElementOfLandauMatrix(int k1,int k2,int i1,int i2);
  
  /*!
    \fn double allLandauMatricesElement(int i_st1,int i_st2,int i,int j)
    \brief provides access to element of landau-matrix
    i_st1:
    i_st2:
    i:
    j:
  */
  std::complex<double> allLandauMatricesElement(int i_st1,int i_st2,int i,int j);

  const yosState *st1,*st2;          

  /*!
    \brief the same as myBasis in Operator
    The really implemented one (yosBasis e.g.)
    And the same after cast (ya like it, Chris, huh? :-)
  */
  yosBasis *endYosBasis; 

  /*!
    \brief 
    How often should the progress
    of setting up the Landau matrix
    should be reported. 
  */
  static const int FREQ_REPORT_OPOP=100; 

  /*!
    \brief number of flux quanta
  */
  int Nm;

  /*!
    \brief number of electrons
  */
  int Ne;

  /*!
    \brief Position to evaluate the Operator at
  */
  double x,y;


  /*!
    \brief When Landau matrices are to be evaluated for the whole basis;
     Numbering: i_st1*dim*Nm*Nm+i_st2*Nm*Nm+i*Nm+j 
  */
  std::complex<double> * m_pCmplxLandauMatrix;

  /*!
    \brief wether allLandauMatrices is already allocated
  */
  bool landau_matrix_allocated;

  /*! \fn whatTypeOfLandauMatrixDoINeed()
      \brief Returns 1 in OneParticleOperator (and DensOperator, jxOp, jyOp). 
      Everytime when getOneElementOfLandauMatrix is rewritten, this
      method should be changed.
  */
  virtual int whatTypeOfLandauMatrixDoINeed();

private:

  
  /*!
    \brief Private methods used by EvalAtXY(double *,int,double,double) and
    matEl_nonzero(). 
  */
  int addToLandauMatrix( int i1, int i2, std::complex<double> *where, std::complex<double> fac );

  double densMatEl(double *coef);

private:
  /*!
    \brief if m_state is valid
  */
  bool m_bStateValid;

  /*!
    \brief state to compute operator for
  */
  eigSt m_stateToCompute;

  /*!
    \brief if this operator needs the wavefunctions to be evaluated
    this flag is needed to dertermine, if the wavefunctions-cahce in
    yosState is to be updated when moving to new position (x,y)
  */
  bool m_bNeedsWF;          


  ///////////////////////////////////////////////////////////////////////
 // needed to share one Landau-Matrix among serveral one-particle-ops //
///////////////////////////////////////////////////////////////////////

  /*!
    \brief if m_pLandauMatrix was already evaluated
  */
  static bool m_bLandauEvaluated;

  /*! 
    \brief what type is the already evaluated Landau matrix; see method
    whatTypeOfLandauMatrixDoINeed();
   */
  static int m_iTypeOfLandauMatrix;

  /*!
    \brief eigSt, for which m_pLandauMatrix is evaluated
  */
  static eigSt m_stateLandauEvaluated;

  /*!
    \brief static pointer to landau matrix, which may be used by other one-particle-operators
  */
  static const std::complex<double> *m_pGlobalLandauMatrix;

protected:
  
  /*!
    \brief Computes s1_nr[i]=2*j[i]+spin[i] for st1. This number may be then used to
    identify a one-electron state uniquely (and also to order the one-electron
    states (using s1_nr[i]<s1_nr[i]). 
  */
  virtual int s1_nr(int i) const;
  /*!
    \brief The same for st2. 
  */
  virtual int s2_nr(int i) const;
};

#endif


