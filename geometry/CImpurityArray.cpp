 #include <geometry/CImpurityArray.h>
#include <LLLlib/LLLlib.h>
#include <utils/CFileParser.h>
#include <utils/CLine.h>

#include <fstream>

#include <geometry/CxPeriodicPosition.hpp>
#include <utils/CDistribution.hpp>
/*
 *  GSL includes
 */
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_histogram2d.h>
 /*
  * System includes
  */
 #include <stdio.h> 
  
CImpurityArray::CImpurityArray():
xSize(1.0f),
ySize(1.0f)
{
	xSize = 1.0f;
	ySize = 1.0f;
}

CImpurityArray::~CImpurityArray()
{
	m_impurities.clear();
}
CImpurityArray::CImpurityArray (const CImpurityArray &rhs)
{
	xSize = rhs.getXSize();
	ySize = rhs.getYSize();
}

CImpurityArray::CImpurityArray(int n_imps, float n_xSize, float n_ySize)
: 	m_impurities(impVector(n_imps)),
	xSize(n_xSize),
	ySize(n_ySize)
{
	
}

int CImpurityArray::addImpurity( const CxImpurity & n_imp)
{
	// check impurity position
	ySize = 1.0f;
	xSize = 1.0f;
	/*if ( (n_imp.getXPosition() > xSize) || (n_imp.getYPosition()> ySize) )
	{
		throw CxErrors("Impurity outside the box");
	}
	*/
	m_impurities.push_back(n_imp);
	return (m_impurities.size());	
	
}

int CImpurityArray::readFromFile(const std::string & fileName)
{
	int retVal = -1;
	CFileParser * myP = CFileParser::getInstance();
	myP->readFile(fileName);
	while (myP->hasNextLine())
	{
		std::vector<float> values =myP->getNextLine().getFloats();
		if (values.size() != 5)
		{
			continue;
		}
		CxImpurity nImp(values[0], values[1], values[2], values[3], values[4]);
		retVal=addImpurity(nImp);
	}	
	retVal = m_impurities.size();
	glLogger.error("File has (%d) impurities", retVal);
	return (retVal);	
}

float CImpurityArray::getPotential(const CxPosition &mPos)const
{	
	float retVal = 0.0f;
	for (unsigned int index = 0; index < m_impurities.size(); index++)
	{
	
		retVal = retVal+m_impurities[index].getPotential(mPos);
		
	}	
	return (retVal);
}

void CImpurityArray::writePositionsToFile(const std::string fileName) const
{
	std::ofstream outStream(fileName.c_str());
	writePositionsToStream(outStream);
	outStream.close();	
	
	
}

void CImpurityArray::writePositionsToStream(std::ostream & outStream) const
{
	outStream <<" # xPos \t \t	yPos \t\t sigmax\t\t sigmay \t\t strenght\n";
	for (unsigned int index = 0; index < m_impurities.size(); ++index) {
		outStream << m_impurities[index].getXPosition()<<"\t";
		outStream << m_impurities[index].getYPosition()<<"\t";
		outStream << m_impurities[index].getSigmaX()<<"\t";
		outStream << m_impurities[index].getSigmaY()<<"\t";
		outStream << m_impurities[index].getStrength()<<"\n";
		
	}
	
}

double CImpurityArray::getExtrema(double &maximum, double &minimum, int gridCount) const
{
	double xEps = 1.0 / gridCount;
	double yEps = 1.0 / gridCount;
	double tempMax = -1e36;
	double tempMin = 1e36;
	for (double xPos = 0.0; xPos < 1.0; xPos = xPos + xEps)
	{
		for (double yPos = 0.0; yPos < 1.0; yPos = yPos + yEps)
		{
			CxPosition temp(xPos,yPos);
			double m_value = getPotential(temp);
			if (m_value > tempMax)
			{
				tempMax = m_value;
			}
			if (m_value < tempMin)
			{
				tempMin = m_value;
			}
				
			
		}
		
	}
	maximum = tempMax;
	minimum = tempMin;
	return (abs(maximum - minimum));	
	
}

void CImpurityArray::writePotentialToFile(std::string m_fileName, int xSteps, int ySteps ) const
{
	std::ofstream outStream(m_fileName.c_str());
	int m_precision = 5;
	if ( (xSteps > 1000) || (ySteps > 1000) )
		{
			m_precision = 7;
		}
	
	outStream.precision(m_precision);
	writePotentialToStream(outStream);
	outStream.close();	
	
}
void CImpurityArray::writePotentialToStream(std::ostream & oStream, int xSteps, int ySteps ) const
{
	if ( (xSteps < 1 ) || (ySteps < 1) )
	{
		return;
	}
	float xEps = 1.0f/xSteps;
	float yEps = 1.0f/ySteps;
	float 
		xPos = 0.0f,
		yPos = 0.0f;
	for (xPos = 0; xPos <= 1.0f; xPos= xPos + xEps) {
		
		for (yPos = 0; yPos <= 1.0f ; yPos= yPos + yEps) {
			CxPosition temp(xPos,yPos);
			oStream << xPos << "  \t  " << yPos << "  \t  "<< getPotential(temp)<<std::endl;
		}
		oStream <<std::endl;
	}
		
		return;
}

std::vector<int> CImpurityArray::createHistogram(int noSteps, int xSteps, int ySteps) const
{
	double 
		max = 0.f,
		min = 0.f,
		width = getExtrema(max,min, xSteps);
	if ( (xSteps < 1 ) || (ySteps < 1) || (noSteps < 1 ))
	{
	//! Throw something here!
		throw CxBadValueError(__FILE__,__LINE__);
	}
	double eps = width / noSteps;
	std::vector<int> retVal(noSteps);
	// loop over all specified points 
	float xEps = 1.0f/xSteps;
	float yEps = 1.0f/ySteps;
	float 
		xPos = 0.0f,
		yPos = 0.0f;
	for (xPos = 0; xPos <= 1.0f; xPos= xPos + xEps) {
		
		for (yPos = 0; yPos <= 1.0f ; yPos= yPos + yEps) {
			CxPosition temp(xPos,yPos);
			double m_potential= getPotential(temp);
			double renorm = m_potential - min; // this is the absolute pot w.r.t. minimum
			unsigned int pos = static_cast<unsigned int> (abs(floor(renorm/ eps)));
			retVal.at(pos) = retVal.at(pos) +1;
		}
		
	}
	return (retVal);
}

CDistribution CImpurityArray::getPairCorrelation(int xSteps, int ySteps) const
{
		
	CDistribution retVal(1000);
	float xEps = 1.0f/xSteps;
	float yEps = 1.0f/ySteps;
	double potentialEps = 1e-06;
	float 
		xPos = 0.0f,
		yPos = 0.0f;
	for (xPos = 0.0f; xPos <= 1.0f; xPos = xPos + xEps)
	{
		for (yPos = 0.0; yPos <=1.0f; yPos = yPos + yEps)
		{
			CxPeriodicPosition aufPunkt(xPos,yPos);
			double m_target= getPotential(aufPunkt);
			/*
			 * Suche Punkte mit dem gleichen Potential (einfach einmal durch)
			 */
			 for (float xSample = 0.0f; xSample <= 1.0f; xSample= xSample+xEps) {
				for (float ySample = 0.0f; ySample <= 1.0f; ySample= ySample+yEps) {
					CxPeriodicPosition testPunkt(xSample,ySample);
					double pot=getPotential(testPunkt);
					if (abs(pot-m_target) < potentialEps)
					{
						float distance = (testPunkt-aufPunkt).vabs();
						retVal.addValue(distance);
						//FOUND
					}
				}
			}
			CxPosition testPunkt;
			
		}
		
		
		
	}
	
	return retVal;
}

bool CImpurityArray::writeHistogramToFile(const char* fileName)const
{
	
	return (true);
}
bool CImpurityArray::write2dHistogramToFile(const char* fileName)const
{
	/*
gsl_histogram2d * h = gsl_histogram2d_alloc (100, 100);
     
       gsl_histogram2d_set_ranges_uniform (h, 0.0, 1.0, 0.0, 1.0);
	for (float aufPunktX = 0.0f; aufPunktX<= 0.0; aufPunktX= aufPunktX+0.01f)
	{
		for (float aufPunktY = 0.0f; aufPunktY<= 0.0; aufPunktY= aufPunktY+0.01f)
		{
		double m_potential = getPotential(CxPosition(aufPunktX, aufPunktY));      
       
			for (float xPos = 0; xPos < 1.0f; xPos = xPos + 0.01) 
			{
				for (float yPos = 0; yPos < 1.0f; yPos = yPos + 0.01) 
				{
					CxPosition temp(xPos,yPos);
					double pot=getPotential(temp);
					if (abs(pot - m_potential) < 1e-04)
					{
				
						gsl_histogram2d_increment (h, xPos, yPos);
					}
				}
			}
		}
	}
    FILE* mFile = fopen(fileName,"w");
    gsl_histogram2d_fprintf (mFile,h,"%f","%f");

    fclose (mFile);
    gsl_histogram2d_free(h);
      */
	return (true);	
}
