#include "LLLRandomLandauMatrix.h"
#include <utils/logger.hpp>
#include <utils/CRandomizer.h>
#include <ERRORS.h>

LLLRandomLandauMatrix::LLLRandomLandauMatrix(
		yosBasis *new_basis,int new_type,double bli_new, 
		double n_corrLength, double n_strength, int n_myImax, double a_new, 
		double b_new):
LLLhamiltonian(new_basis, new_type, bli_new, a_new, b_new),
m_uOfk(std::vector<double>()),
m_iMax(abs(n_myImax)),
m_corrLength(n_corrLength),
m_strength(n_strength),
m_sigma(getSecondExponent()),
m_prefactor(0.0), 
m_coefficents_dim(0),
m_coefficents(0)
{
	m_coefficents_dim = (new_basis->getNm()*2+m_iMax)*(new_basis->getNm());
	m_coefficents = new std::complex<double>[m_coefficents_dim];
	//m_coefficents = new std::complex<double> [dim *dim]);
	calcPreFactors();
	populateMatrix();
}



LLLRandomLandauMatrix::~LLLRandomLandauMatrix()
{
	delete [] m_coefficents;
	m_coefficents = 0;
}


/*! \fn double LLLRandomLandauMatrix::getSecondExponent()const
 \brief calculates \sigma in Huckestein  
 * 
 */
double LLLRandomLandauMatrix::getSecondExponent()const
{
	double retVal = 0.0;
	
	for (int m_i = 0; m_i < m_iMax; ++m_i) 
	{
		retVal = retVal + getExpOfI(2*m_i);
		
		
	}	
	
	return (retVal);
}

double LLLRandomLandauMatrix::getExpOfI(const int index )const
{
	int Nm = static_cast<yosBasis*> (my_basis)->getNm();
	double prefac = -M_PI/ Nm;
	prefac = prefac / m_corrLength / m_corrLength;
	prefac = prefac * index * index;
	prefac = exp(prefac);	
	return (prefac);
	
}
/*!
 \brief calculates everything before the u dependent sum
 */
void LLLRandomLandauMatrix::calcPreFactors(void)
{
	double a_to_b = get_Lx()/get_Ly();
	double preFac = M_PI*2;
	double rtNm = sqrt((double)static_cast<yosBasis*> (my_basis)->getNm());
	preFac = preFac*rtNm*m_sigma;
	preFac = preFac *a_to_b;
	preFac = sqrt(preFac);
	preFac = m_strength/preFac;
	for (int deltak = 0; deltak < Nm; ++deltak) 
	{
		int deltak2 = deltak*deltak;
		double exponent = deltak2*m_corrLength*m_corrLength*0.25;
		exponent = exp(exponent);
		deltak_Prefacs.push_back(exponent*preFac);
		
	}

	
} 



void LLLRandomLandauMatrix::setCorrelationLength(double new_corrLength)
{
	m_corrLength= new_corrLength;	
	
}

int LLLRandomLandauMatrix::matEl_perturb_nonzero(int iIndex,int jIndex ,std::complex<double> *retVal)
{
	
	int 
		rIndex = 0, 
		lIndex = 0,
		sign = 1;
		glLogger.debug("Entering LLLRandomLandauMatrix::matEl_perturb_nonzero with (%i), (%i)", iIndex, jIndex);
	st1=(*LLL_basis)[iIndex];
	st2=(*LLL_basis)[jIndex];
	*retVal = std::complex<double>(0.0,0.0);
	if (iIndex == jIndex)
	{
		// diagonal element	sum over all elements
		for (int eIndex = 0; eIndex < st1->getNe(); ++eIndex) 
		{
			std::complex<double> tempVal =  getDiagonalElements(st1->j(eIndex));
			*retVal = *retVal + tempVal;
		}
	} 
	else 
	{
		bool control = 	st1->getSingleDiffOccNumber(*st2,rIndex, lIndex, sign);
		if (control)
		{
			*retVal = getOffDiagonalElements(st1->j(rIndex),st2->j(lIndex))* std::complex<double>(sign, 0);
		}
		else 
		{
			glLogger.debug("null returned");
		}
	}

	if (fabs(retVal->imag()) < 1e-20) 
	{
		std::complex<double> temp = std::complex<double>(retVal->real(), 0.0);
		retVal = &temp;
	}
	glLogger.debug("matelPerturb returns (%f)", retVal->real());
return (1);	
	
}
/*!
 * \brief calculates diagonal elements  
 */
std::complex <double> LLLRandomLandauMatrix::getDiagonalElements(int index) const
{
	std::complex <double> retVal;
	glLogger.debug("LLLRandomLandauMatrix::getDiagonalElements (%d)", index);
	retVal = deltak_Prefacs.at(0);
	for (int iIndex = 0; iIndex < m_iMax; ++iIndex) {
		retVal = retVal + get_U_ofK(2*index+iIndex, 0)*getExpOfI(iIndex);
	}
	std::complex <double> tretVal(retVal.real(), 0.0);
	return (tretVal);	
	
}

std::complex<double> LLLRandomLandauMatrix::getOffDiagonalElements(int iIndex, int jIndex) const
{
	std::complex <double> retVal, temp;
	int deltak = abs(iIndex -jIndex);
	retVal =deltak_Prefacs.at(deltak);
	for (int m_i = 0; m_i < m_iMax; ++m_i) 
	{
		temp = temp+get_U_ofK(2*iIndex+m_i,deltak)*getExpOfI(m_i);
	}	
		retVal = std::complex<double>(temp.real(),0.0);
	return (retVal);	
	
	
}


bool LLLRandomLandauMatrix::populateMatrix()
{

	if (!m_coefficents) 
	{
		return (false);
	}
	CRandomizer *myRan = CRandomizer::getInstance();
	myRan->randomize(198912);
	// calculate the dimension of the array
	for (int m_index = 0; m_index < m_coefficents_dim; ++m_index) {
		m_coefficents[m_index] = std::complex<double>(myRan->dRand(-m_strength,m_strength), 0.0);
		glLogger.debug("Random Matrix [%d] = (%f)", m_index, m_coefficents[m_index].real());
	}
	
	// loop  to fill the matrix 
		
	return (true);
}
std::complex<double> LLLRandomLandauMatrix::get_U_ofK(int row, int col) const
{
	if (row <0 || col < 0)
	{
		throw CxIndextoSmallError(__FILE__, __LINE__); 
	}
	std::complex<double> retVal;	
	glLogger.info("Random Matrix element (%d),(%d)", row, col);	
	
	int tempIndex = row*(LLL_basis->getNm()* 2 * m_iMax) +col;
	try {
		retVal= m_coefficents[tempIndex];
		glLogger.info("Random Matrix element (%d),(%d) = (%f)", row, col, retVal.real());	
	}
	catch (std::exception e)
	{
		glLogger.error("Trying to read out of bounds (%d),(%d), %s", row, col, e.what());
		throw CxIndextoSmallError(__FILE__, __LINE__); 
	}
	catch (...)
	{
		glLogger.error("caught something else in  LLLRandomLandauMatrix::get_U_ofK");	
		throw CxErrors(__FILE__, __LINE__);
	}
	return (retVal);
}
