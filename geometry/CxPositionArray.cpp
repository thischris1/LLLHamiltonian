#include <geometry/CxPositionArray.h>
#include <math.h>
CxPositionArray::CxPositionArray():
  m_positions(std::vector<CxPosition>())
{

}

CxPositionArray::~CxPositionArray()
{
  m_positions.clear();
}


/*! 
  \brief Search if position is already in the array.
  \return the index of the position, -1 if not found
  \param nPos the position to look up
*/


int CxPositionArray::findPosition(const CxPosition &nPos) const
{
  int retVal = -1;
 
  for (unsigned int index = 0; index < m_positions.size(); index++)
    
    {
      CxPosition m_pos = m_positions[index];
      if (nPos.compare(m_pos))
	{
	  retVal = index;
	  break;
	}

    }
  return (retVal);
}

bool CxPositionArray::addPosition( const CxPosition &n_Pos)
{
  if (findPosition(n_Pos) ==-1)
    {
      return (false);
    }
  m_positions.push_back(n_Pos);
  return (true);

}

std::vector<CxPosition> CxPositionArray::getAllAroundPos(const CxPosition& mittelPunkt, float distance)const
{
	std::vector<CxPosition> retVal;
	float m_distance = fabs(distance);
	// loop over all positions and calculate distance
	
	for (unsigned int index = 0; index < m_positions.size(); index++)
	{
	 	CxPosition m_pos = m_positions[index];
		float tempDistance = (m_pos-mittelPunkt).vabs();
		if (tempDistance < m_distance)
		{
			// belongs to
			retVal.push_back(m_pos);
		}
		
		
	}	
	return (retVal);
	
}
