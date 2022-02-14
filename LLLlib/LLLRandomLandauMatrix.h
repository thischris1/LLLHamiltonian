#ifndef LLLRANDOMLANDAUMATRIX_H_
#define LLLRANDOMLANDAUMATRIX_H_

#include "LLLhamiltonian.hpp"
#include <vector>


class LLLRandomLandauMatrix : public LLLhamiltonian
{
public:
	LLLRandomLandauMatrix(yosBasis *new_basis,int new_type,double bli_new,
		double n_corrLength, double n_strength, 
		int n_myImax=20, double a_new=1.0, double b_new=1.0);
	virtual ~LLLRandomLandauMatrix();
	
	void setIMax(int n_imax);
	
	int getIMax(void)const {return m_iMax;};
	
	double getCorrelationLength(void)const {return (m_corrLength);}; 
	
	void setCorrelationLength(double new_corrLength);
	
	double getStrength(void) const {return (m_strength);};
	void setStrength(double n_strength){m_strength= n_strength;};
	virtual int matEl_perturb_nonzero(int,int,std::complex<double> *);
	
	double getSigma(void) const {return m_sigma;};
	
	std::vector<double> getDeltak_prefacs(void ) const { return (deltak_Prefacs);};

protected:
//! helper function
	std::complex<double> getDiagonalElements(int index) const;

//! helper function
	std::complex<double> getOffDiagonalElements(int iIndex, int jIndex) const;
	
	
	
	
private:
	// helper functions
	//! calculate sum over exponentials \approx k_{1} - k_{2}
	double getFirstExp(int delta_k) const;
	//! calc. the second exponent (factor in sum and norm. from huckesteins formual)
		
	double getSecondExponent()const;
	
	void calcPreFactors();
	//! Helper for the sum of exp 
	double getExpOfI(const int index) const; 	
	//! Populates the random number matrix
	bool populateMatrix();	//! 
//! Accessor to internal array 
	std::complex<double> get_U_ofK(int row, int col) const;
	
	// data members
	std::vector<double> m_uOfk; //! list of randomly generated values
	int m_iMax;  //! index for summation
	double m_corrLength; //! the factor \beta in Huckestein
	double m_strength; //! V_{0} in Huckestein
	double m_sigma ; //! Prefactor \Sigma from Huckestein's paper
	double m_prefactor; //! \frac{V_{0}}{{\sqrt(2) L_{y}* \Sigma)^{\frac{1}{2}}
	int m_coefficents_dim; //! dimension of array m_coefficents
	std::vector<double> deltak_Prefacs;
	std::complex<double> *   m_coefficents; //! u array from Huckestein
};

#endif /*LLLRANDOMLANDAUMATRIX_H_*/
