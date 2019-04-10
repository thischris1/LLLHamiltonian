#include "CGaussianImpurityArray.h"
#include <utils/CRandomizer.h>
#ifdef GSL
#include <utils/CRandomizerGsl.h>
#endif 
#include <geometry/CxPosition.hpp>
#include <utils/logger.hpp>
#include <utils/CFileParser.h>
#include <utils/CLine.h>

#include <fstream>


CGaussianImpurityArray::CGaussianImpurityArray(int  noOfImpurities, float  sigma , float n_maxStrength, float n_minStrength, float n_xSize, float n_ySize):
CImpurityArray(0, n_xSize, n_ySize),
m_strengthLimitMax(n_maxStrength),
m_strengthLimitMin(n_minStrength),
m_StrengthSum(0.0)
{
	glLogger.debug("Creating a gaussian impurtiy array with size (%d), sigma(%f)", noOfImpurities, sigma);
	std::vector <CxPosition> m_pos(noOfImpurities);

	#ifdef GSL
	if (CRandomizerGsl::getInstance()->randomPositions(m_pos) !=  noOfImpurities)
	{
		return;
	}
	#else
	if (CRandomizer::getInstance()->randomPositions(m_pos) !=  noOfImpurities)
	{
		return;
	}
	#endif
	 // add too array
	 
	 int sign = -1;
	for (unsigned int index = 0; index < (unsigned int)noOfImpurities; ++index) 
	{
		#ifdef GSL
		float tempStrength = CRandomizerGsl::getInstance()->fRand(m_strengthLimitMin, m_strengthLimitMax);
		#else
			float tempStrength = CRandomizer::getInstance()->fRand(m_strengthLimitMin, m_strengthLimitMax);
			#endif
//		CxImpurity tempImp(m_pos[index], sigma, sigma, strength*sign);
		CxImpurity tempImp(m_pos[index], sigma, sigma, tempStrength);	
		sign= sign*-1;
		addImpurity(tempImp);
		m_StrengthSum+= tempStrength;
	}
	 
	
}


/*! \brief constructor reading from file
 * 
 */

CGaussianImpurityArray::CGaussianImpurityArray(const std::string &fileName):
m_strengthLimitMax(1.0),
m_strengthLimitMin(0.0),
m_StrengthSum(0.0)
{

	CFileParser *m_parser = CFileParser::getInstance();
	m_parser->readFile(fileName);
	int noOfImpurities = m_parser->getNextLine().getIntegers().at(0);
	float sigma = m_parser->getNextLine().getFloats()[0]; 
	m_strengthLimitMax = m_parser->getNextLine().getFloats()[0];	
	m_strengthLimitMin = m_parser->getNextLine().getFloats()[0];


	std::vector <CxPosition> m_pos(noOfImpurities);

	if (CRandomizer::getInstance()->randomPositions(m_pos) !=  noOfImpurities)
	{
		return;
	}
		
	 // add too array
	for (unsigned int index = 0; index < (unsigned int)noOfImpurities; ++index) 
	{
		#ifdef GSL
		CxImpurity tempImp(m_pos[index], sigma, sigma, CRandomizerGsl::getInstance()->fRand(m_strengthLimitMin, m_strengthLimitMax));
		#else
		CxImpurity tempImp(m_pos[index], sigma, sigma, CRandomizer::getInstance()->fRand(m_strengthLimitMin, m_strengthLimitMax));
		#endif
		addImpurity(tempImp);
	}
	 
	
}

CGaussianImpurityArray::~CGaussianImpurityArray()
{
}

void CGaussianImpurityArray::writePositionsToFile(const std::string fileName) const
{
	std::ofstream outStream(fileName.c_str());
	outStream <<" # random Array \n";
	CImpurityArray::writePositionsToStream(outStream);
	
	outStream.close();	
	
	
}

bool CGaussianImpurityArray::setStrengthLimits(const double & n_strengthMax, double & n_strengthMin)
{
	if (n_strengthMax < n_strengthMin)
		{
			return (false);
		}
	m_strengthLimitMax = n_strengthMax;
	m_strengthLimitMin = n_strengthMin;
	// populate array with new psoitions
	
	
	return (true);	
} 

