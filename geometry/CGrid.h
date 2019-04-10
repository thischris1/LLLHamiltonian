#ifndef _CGRID_H
#define _CGRID_H
#include <geometry/CxPosition.hpp>

class CGrid {
  public:
    CGrid();

    CGrid(double xSize, double ySize, int xSteps, int ySteps);

    virtual ~CGrid();

    // Accessors
    inline double getXLength() const {return (m_xLength);};

    inline double getYLength() const {return (m_yLength);};

    inline double getYMeshSize() const {return (m_yMeshSize);};

    inline double getXMeshSize() const {return (m_xMeshSize);};

    inline int getxSteps() const {return (m_xSteps);};

    inline int getySteps() const {return (m_ySteps);};
	

  protected:
    int getNearestXIndex(const double & xPoint) const;

    int getNearestYIndex(const double & yPoint) const;

    double getNearestXPoint(const double & xPoint) const;

    double getNearestYPoint(const double & yPoint) const;

//! We assume the value to be stored in one long array, row-wise
	int mapIndex (int xPoint, int yPoint) const;
	
	int mapDoubleToIndex (double xPoint, double yPoint) const;
	int positionToIndex (const CxPosition &n_Pos) const;
	CxPosition indexToPosition(int rowIndex, int colIndex) const;

	int meshIndexToCoordinate(int x, int y, double& new_x, double & new_y)const;
  private:
    double m_xLength;

    double m_yLength;

    int m_xSteps;

    int m_ySteps;

    double m_xMeshSize;

    double m_yMeshSize;

/*
 *  organization of mesh is as follows:
 *  zeilenweise von 0 ... xSteps
 *  0 _mySteps
*/
};
#endif
