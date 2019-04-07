#include "CDistanceDistributionHR.h"
#include <LLLlib/LLLlib.h>
#include <math.h>
#include <fstream>

CDistanceDistributionHR::CDistanceDistributionHR():
m_minVal(0.0f),
m_maxVal(1.0f),
m_startHighRes(0.0f),
m_endHighRes(0.0f),
m_lowResolution(1000),
m_highResolution(1000)
{
	// Generates a uniform distribution
}
/*!
 \param n_minVal
 \param n_maxVal
 \param n_lowRes
 \param n_highResStart
 \param n_highRes
 \param n_highResEnd
 \brief Ctor for a Distribution
 * 
 * 
 */

CDistanceDistributionHR::CDistanceDistributionHR(float n_minVal, float n_maxVal, unsigned int n_lowRes, unsigned int n_highRes, float n_highResStart, float n_highResEnd):
CDistribution(n_lowRes+n_highRes),
m_minVal(n_minVal),
m_maxVal(n_maxVal),
m_startHighRes(n_highResStart),
m_endHighRes(n_highResEnd),
m_lowResolution(n_lowRes),
m_highResolution(n_highRes)
{
	// check sanity
	if ((n_minVal > n_maxVal) || fabs (n_maxVal - n_minVal) < 1e-04)
	{
			throw CxBadValueError(__FILE__,__LINE__, "bad Limits of distribution");
	}
	if ( (n_highResEnd < n_highResStart ) || fabs (n_highResEnd - n_highResStart) < 1e-05)
	{
		
		throw CxBadValueError(__FILE__,__LINE__, "bad Limits of high resolution");
	}
	if ( (n_highResStart < n_minVal ) || (n_highResEnd > n_maxVal))
	{	
		  throw CxBadValueError(__FILE__,__LINE__, "bad Limits of high resolution (too large)");
	}	
	if ( (n_lowRes < 2) ||   (n_highRes) < 2 )
	{
		throw CxBadValueError(__FILE__,__LINE__, "bad Limits of high resolution (too large)");
	}
	// calculate intervall sizes now
	float lengthOfHighIntervall =  fabs (n_highResEnd - n_highResStart);
	m_highResInterVal = lengthOfHighIntervall / m_highResolution;
	// low length = total length - highRes Length
	float lengthOfLowIntervall = n_maxVal - n_minVal - lengthOfHighIntervall;
	m_lowResInterVal =  lengthOfLowIntervall / m_lowResolution;
	m_Values = std::vector<int>(n_lowRes+n_highRes);
	// find start index for high resolution
	if (fabs(m_startHighRes - m_minVal) < m_highResInterVal)
	{
		m_startIndexHighResolution = 0;	
		m_endIndexHighResolution = m_highResolution;
		
	}
	else
	{
		// calculate count lower Part / total low resolution
		float tempLowPart = (m_startHighRes - m_minVal) /  lengthOfLowIntervall;
		float tempUpperPart = (m_maxVal - m_endHighRes) / lengthOfLowIntervall;
		if ( fabs(tempLowPart + tempUpperPart -1.0f) > 1e-04)
		{
			throw   CxBadValueError(__FILE__,__LINE__, "Could not calculate limits");
		}
		m_startIndexHighResolution = static_cast<size_t>(floor(tempLowPart * m_lowResolution)) + 1;
		m_endIndexHighResolution = m_startIndexHighResolution + m_highResolution;
	}
	
}

CDistanceDistributionHR::~CDistanceDistributionHR()
{
	m_Values.clear();
}


CDistanceDistributionHR::CDistanceDistributionHR (const CDistanceDistributionHR &rhs ):CDistribution(rhs)
{
	if (&rhs	 != this)
	{
		m_minVal = rhs.getMinimumValue();
		m_maxVal = rhs.getMaximumValue();
		m_startHighRes = rhs.getStartHighResolution();
		m_endHighRes = rhs.getEndHighResolution();
		m_Values = rhs.getHistogram();		
	}

	
	
}

bool CDistanceDistributionHR::addValue(const double n_value)
{
	if (!isValidValue(n_value) )
	{
		return (false);
	}
	
	if (isInHighResolutionPart(n_value))
	{
		// Calculate index+
		double tempVal = n_value - m_startHighRes;
		unsigned int tempIndex = static_cast<unsigned int>(fabs(round (tempVal / m_highResInterVal)));
		unsigned int fullIndex = m_startIndexHighResolution + tempIndex;
		m_Values.at(fullIndex) = m_Values.at(fullIndex)+1;
		
	}
	else 
	{
		
		/* decide whether it is in the intervall before or after the 
		 * highResoultion interval
		 */
		 if ( (n_value < m_startHighRes) )
		 {
		 	// is valid therefore not smaller
		 }
		 else 
		 {
		 	//after highres
		 	//unsigned int lowCount = m_lowResolution;
		 	// calculate distance to end of high Res
		 	double tempDistance = n_value - m_endHighRes;
		 	// calculate pocket count
		 	unsigned int tempIndex = static_cast<unsigned int>(fabs(round(tempDistance / m_lowResInterVal)));
		 	unsigned int finalIndex = m_startIndexHighResolution + m_highResolution + tempIndex;
		 	m_Values.at(finalIndex) = m_Values.at(finalIndex)+1;
		 }  
		
		
	}
	return (true);	
}


bool CDistanceDistributionHR::addValue(const float n_val)
{
	return (addValue(static_cast<double>(n_val) ) );
}



bool CDistanceDistributionHR::isValidValue(const double n_val)
{
	if ( (n_val < m_minVal) << (n_val > m_maxVal) )
	{
		return (false);
	}
	return (true); 	
}
/*! 
 * \brief convenience method to distinguish  low res from high res positions
 * 
 */
bool CDistanceDistributionHR::isInHighResolutionPart(const double n_val)
{
		if (!isValidValue(n_val))
		{
			return (false);
		}
		
		if ( (m_startHighRes < n_val) && (m_endHighRes >  n_val) )
		{
			return (true);
		} 
		else 
		{
			return (false);
		}
}

bool CDistanceDistributionHR::writeToFile(const std::string & n_fileName) const
{
	std::ofstream outStream(n_fileName.c_str());
	if (outStream.good())
	{
		bool retVal = writeToStream(outStream);
		outStream.close();
		return (retVal);
	}
	else
	{
		return (false);
	}	
}

bool CDistanceDistributionHR::writeToStream(std::ostream &outStream) const
{
	if (!outStream.good())
	{
			
		return (false);
	}
	// write header
	outStream << " # Distribution with different resolutions \n";
	outStream << " # high Resolution starts at " << m_startHighRes << " at index (" << m_startIndexHighResolution;
	outStream << ") ends at " << m_endHighRes<< " index ("<< m_endIndexHighResolution << ")" <<std::endl;
	outStream << "# Minimum / Start value = " << m_minVal << " end at " << m_maxVal << std::endl;
	outStream << "# low Resolution = " << m_lowResolution  << " High resolution = "<< m_highResolution <<std::endl; 
	outStream << "# index  \t radius \t count 1 = High resolution 0 = lowResolution\n";
	outStream << std::endl;
	float position = 0.0;
	for (unsigned int index = 0; index < m_Values.size(); ++index) {
		// Calculate distance
		bool inHigh = false;
		if (index == 0)
		{
			// do not count up
		}
		else 
		{
			if (index < m_startIndexHighResolution)
			{
				// in low resolution part 
				position = position + m_lowResInterVal;
				inHigh = false;
			}
			else if (index > m_endIndexHighResolution)
			{
				// after high Resolution part
				position = position + m_lowResInterVal;
				inHigh = false;
			}
			else 
			{ 
				// in high resolution part 
				position = position + m_highResInterVal;
				inHigh= true;
			}
	// 	
		}
		outStream << index << "\t " << position <<"\t "  << m_Values.at(index) << "\t "<< inHigh<< std::endl;
	}
	return (true);
}

