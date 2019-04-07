#include <utils/C2dDistribution.hpp>
#include <math.h>

C2dDistribution::C2dDistribution():
  N_STEPS_X(1000),
  N_STEPS_Y(1000),
  m_Values(C2DVector<int>(N_STEPS_X,N_STEPS_Y))
{

}
C2dDistribution::C2dDistribution(const int & x_size, const int ySize):
  N_STEPS_X(x_size),
  N_STEPS_Y(ySize),
  m_Values(C2DVector<int>(N_STEPS_X,N_STEPS_Y))
{

}


C2dDistribution::~C2dDistribution()
{
}

bool C2dDistribution::writeToStream(std::ostream &m_stream) const
{
  if (!m_stream.good())
    {
      return (false);
    }
  for (int xIndex=0; xIndex < (int)N_STEPS_X; xIndex++)
    {
      for ( int yIndex=0; yIndex < (int)N_STEPS_Y; yIndex++)
	{
	  int value =  m_Values.GetAt(xIndex, yIndex);
	  m_stream << xIndex<< "  "<<yIndex << "  " << value <<std::endl;
	}
    }

  return (true);
}


bool C2dDistribution::addValue (const CxPeriodicPosition &newPos)
{
  float xSize = newPos.get_xCellSize();

  int newX = (int)rint( newPos.getXPosition()/xSize * N_STEPS_X);
  int newY = (int)rint( newPos.getYPosition()/newPos.get_yCellSize() * N_STEPS_Y);
  m_Values.IncAt(newX, newY, 1);
  return (true);



}

C2DVector<int>  C2dDistribution::getValues(void) const 
{ 
  return (m_Values);
}



