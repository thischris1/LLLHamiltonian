#ifndef CINTEGERGRID_H_
#define CINTEGERGRID_H_

#include "CScalarGrid.h"
#include <vector>
#include <iostream>

typedef std::vector<int> intRowType;
typedef std::vector<int>::iterator intRowTypeIt;


class CIntegerGrid : public CScalarGrid
{
public:
	
	CIntegerGrid(float xSize, float ySize , int xMeshCount, int yMeshCount);
	virtual ~CIntegerGrid();
	double getValue(const CxPosition &n_pos)const;

	bool setValue (const CxPosition &n_pos, double newValue);	
	bool getExtrema(double & max, double &min) const;
    virtual bool dumpFullContent(void);
protected:
	CIntegerGrid();
	bool initialize();
private:
	//std::vector<std::vector<int > > m_values;
	intRowType m_values;
};

#endif /*CINTEGERGRID_H_*/
