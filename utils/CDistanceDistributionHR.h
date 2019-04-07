#ifndef CDISTANCEDISTRIBUTIONHR_H_
#define CDISTANCEDISTRIBUTIONHR_H_

#include <utils/CDistribution.hpp>

/*!
 \class CDistanceDistributionHR
 \brief class for different resolution of histograms
 The total number of pockets is given by the sum of highres and lowres in the ctor. The start 
 and the end of the high resolution part can be specified. The intervalls are evenly spaced 
 in the high and low res aprt 
 
 */

class CDistanceDistributionHR : public CDistribution
{
public:
	CDistanceDistributionHR(float n_minVal, float n_maxVal, unsigned n_lowRes, unsigned int n_highRes, float n_highResStart, float n_highResEnd);
	virtual ~CDistanceDistributionHR();
	
	CDistanceDistributionHR (const CDistanceDistributionHR &);
	CDistanceDistributionHR operator = (const CDistanceDistributionHR &);
	
	CDistribution normalize(const int & newCount );
	virtual bool addValue(const float);
  	virtual bool addValue (const double);
  	virtual bool writeToStream(std::ostream &) const;
  	virtual bool writeToFile(const std::string & n_fileName) const;
  	
  	float getMinimumValue(void) const {return m_minVal;};
  	float getMaximumValue(void) const {return m_maxVal;};
  	float getStartHighResolution(void) const {return m_startHighRes;};
  	float getEndHighResolution(void) const {return m_endHighRes;};
  	unsigned int getLowResolution(void) const {return m_lowResolution;};
  	unsigned int getHighResolution(void) const {return m_highResolution;};
  	std::vector<int> getHistogram(void) const {return m_Values;};
  	
  	
protected:
	bool isInHighResolutionPart( const float );
	bool isInHighResolutionPart( const double );
	bool isValidValue ( const double );

private:
// no implementation 
	CDistanceDistributionHR();
	//! member variables next 
	float m_minVal;
	float m_maxVal;
	float m_startHighRes; //! The start value (distance) for the higher resolution)
	float m_endHighRes; //! The end value (distance) for the higher resolution
	unsigned int m_lowResolution; //! number of pockets in the low resolution part 
	unsigned int m_highResolution; //! number of pockets in the higher resolution part
	std::vector<int> m_Values;
	float m_highResInterVal;
	float m_lowResInterVal;
	size_t m_startIndexHighResolution;
	size_t m_endIndexHighResolution;
};

#endif /*CDISTANCEDISTRIBUTIONHR_H_*/
