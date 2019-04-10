#ifndef CGAUSSIANIMPURITYARRAY_H_
#define CGAUSSIANIMPURITYARRAY_H_

#include <geometry/CImpurityArray.h>
#include <string>


class CGaussianImpurityArray : public CImpurityArray
{
public:
	CGaussianImpurityArray(int  noOfImpurities, float  sigma, float n_maxStrength = 0.5f, float n_minStrength = -0.5f, float n_xSize = 1.0f, float n_ySize = 1.0f);
	CGaussianImpurityArray(const std::string &fileName); 
	CGaussianImpurityArray operator = (CGaussianImpurityArray & rhs);
	virtual ~CGaussianImpurityArray();
	double getStrengthMax(void) const { return (m_strengthLimitMax);}
	double getStrengthMin(void) const { return (m_strengthLimitMin);}
	double getPotentialAverage(void) const { return (m_StrengthSum);}
	bool setStrengthLimits(const double & n_strengthMax, double & n_strengthMin); 
   void writePositionsToFile(const std::string fileName) const;
   	bool generateHistogram(std::vector<int> m_histogram)const ;
	bool writeHistogramToFile( const std::string &fileName) const;  	
 
   
private:
	CGaussianImpurityArray();
	
	// member variables
	double m_strengthLimitMax;
	double m_strengthLimitMin;
	double m_StrengthSum; // gives sum over all potentials (absolut values)
};

#endif /*CGAUSSIANIMPURITYARRAY_H_*/
