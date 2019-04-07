#ifndef LLLHAMGAUSSIANARRAY_H_
#define LLLHAMGAUSSIANARRAY_H_

#include <LLLhamiltonian.hpp>
#include <yosBasis.hpp>
//#include <geometry/CxImpurity.hpp>
#include <vector>
#include <Types.hpp>
//#include <numerics/SquareMatrix.h>

class LLLhamGaussianArray : public LLLhamiltonian
{
public:
	LLLhamGaussianArray(yosBasis* new_basis,int new_type,double bli_new, double a_new,double b_new, double offset);
	LLLhamGaussianArray(yosBasis* new_basis,int new_type,double bli_new, double a_new,double b_new, double offset, impVector n_imps);
	virtual ~LLLhamGaussianArray();
	virtual int matEl_perturb_nonzero(int,int,CMDCOMPLEX *);
	impVector  getImpurities()const { return (m_impurities);}
	//! add an impurity
	bool addImpurity(const CxImpurity & n_imp);
	//! clears all impurities
	bool clearImpurities(void);
	// write impurities to logger
	bool logImpurities(void) const; 
private:
// only declarations!
	LLLhamGaussianArray();
	LLLhamGaussianArray( const LLLhamGaussianArray &rhs);
	LLLhamGaussianArray operator = (const LLLhamGaussianArray &rhs);
// private helpers
	CMDCOMPLEX getDiagonalElement(const int & index) const;
	//CMDCOMPLEX getOffDiagonalElement(const int& row, const int & column, double alpha1, double alpha2) const;
	CMDCOMPLEX getSPMatrixElement(const int& row, const int & column, double alpha1, double alpha2) const;
	// just to calculate ... 
	CMDCOMPLEX getSPAMAtrixElement(const int& row, const int & column, double alpha1, double alpha2) const;
	MDOUBLE iPlusj (int  i_plus_j, double  alpha2, double &xEvenPass, double &yEvenPass) const;
	CMDCOMPLEX iMinusj(CMDCOMPLEX & MEYEvenPass, CMDCOMPLEX & MEYOddPass, double alpha1, int i_minus_j) const;
// members
	std::vector<CxImpurity> m_impurities;
	MDOUBLE m_offSet; //! to provide an offset energy for diagonalization purposese
	int xKmax;
	int yKmax;
	bool preCalculateMatEls();
	
	CMDCOMPLEX **m_tempArray;  
 

};

#endif /*LLLHAMGAUSSIANARRAY_H_*/
