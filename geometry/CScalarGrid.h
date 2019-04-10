#ifndef CSCALARGRID_H_
#define CSCALARGRID_H_
#include <vector>
#include "CGrid.h"
#include <geometry/CxPosition.hpp>
#include <utils/CDistanceDistribution.h>
#include <utils/IValueProvider.h>
class CScalarGrid : public CGrid
{
public:
	
	CScalarGrid(float xSize, float ySize , int xMeshCount, int yMeshCount);
	virtual ~CScalarGrid();
	/*
	 * Get functions here
	 */
	virtual double getValue(const  CxPosition &n_pos)const = 0;
    double getValueAtPosition(const CxPosition &n_pos)const;
    double getValueAtIndex(const int &, const int &) const;
   /*
    * 
    * Set functions
    */
    CxPosition findMaximumPosition(void)const;
    CxPosition findMinimumPosition(void) const;
    virtual double getMaximum(void)const;
    virtual double getMinimum(void)const;
    virtual bool getExtrema(double & max, double &min)const = 0;
    virtual bool setValue (const CxPosition &n_pos, double newValue) = 0;
    
    // fill the grid with values from Valueprovider
   	bool fillGrid(IValueProvider *m_provider);
   
    CDistanceDistribution  getCorrelation(const int Steps, double eps = 1e-05);
	bool dumpContent(void) const;
	virtual bool dumpFullContent(void) = 0;

protected:
	CScalarGrid();
	bool initialize();
private:
	
protected:
};

#endif /*CSCALARGRID_H_*/
