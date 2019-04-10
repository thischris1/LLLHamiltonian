#ifndef CDOUBLEGRID_H_
#define CDOUBLEGRID_H_
#include <vector>
#include "CScalarGrid.h"

class CDoubleGrid : public CScalarGrid
{
public:
	CDoubleGrid(float xSize, float ySize , int xMeshCount, int yMeshCount);
	virtual ~CDoubleGrid();
	double getValue(const CxPosition &n_pos)const;
	bool setValue (const CxPosition &n_pos, double newValue);	
	bool getExtrema(double & max, double &min) const;
    virtual bool dumpFullContent(void);
protected:
	CDoubleGrid();
	bool initialize(double newVal=0.0);
private:
	std::vector<double > m_values;
};

#endif /*CDOUBLEGRID_H_*/
