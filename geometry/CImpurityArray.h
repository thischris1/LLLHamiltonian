#ifndef CIMPURITYARRAY_H_
#define CIMPURITYARRAY_H_
#include <geometry/CxImpurity.hpp>
#include <string>
#include <LLLlib/LLLlib.h>
#include <LLLlib/Types.hpp>
#include <geometry/CxPosition.hpp>
#include <utils/CDistribution.hpp> 
#include <geometry/CxPositionArray.h>
/*!
 \class CImpurityArray
 */
class CImpurityArray: public CxPositionArray
{
public:
	CImpurityArray();
	//! create n_size random impurities
	CImpurityArray(int n_imps, float n_xSize, float n_ySize);
	virtual ~CImpurityArray();
	// Copy ctor
	CImpurityArray (const CImpurityArray &rhs) ;
	CImpurityArray (const std::string & fileName); 
	impVector getImpurities(void) const { return (m_impurities);};
	float getXSize()const { return xSize;};
	float getYSize() const { return ySize;};
	void setXSize(float n_size) { xSize = n_size;};
	void setYSize(float n_size) { ySize = n_size;};
	//! get the potential at a position (superposition of impurities)
	float getPotential (const CxPosition & mPos) const;
	int addImpurity (const CxImpurity & n_imp);
	int getSize() const { return (m_impurities.size());};
	virtual int readFromFile (const std::string & fileName);
	/* 
	 * Evaluation methods
	 */
	virtual 	void writePositionsToFile(const std::string fileName) const;
	virtual void writePositionsToStream( std::ostream & oStream) const;
	void writePotentialToFile(std::string m_fileName, int xSteps = 100, int ySteps =100 ) const;
	void writePotentialToStream(std::ostream & oStream, int xSteps = 100, int ySteps =100 ) const;
	double getExtrema(double &maximum, double &minimum, int gridCount = 100)const;
	std::vector<int> createHistogram(int noSteps, int xSteps = 100, int ySteps = 100)const;	
	CDistribution getPairCorrelation(int xSteps = 100, int ySteps = 100) const;
	bool writeHistogramToFile(const char* fileName)const;
	bool write2dHistogramToFile(const char* fileName)const;
private:
	impVector m_impurities;
	float xSize;
	float ySize;
	
};

#endif /*CIMPURITYARRAY_H_*/
