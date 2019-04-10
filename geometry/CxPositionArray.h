#ifndef CXPOSITIONARRAY_HPP
#define CXPOSITIONARRAY_HPP
#include <vector>
#include <geometry/CxPosition.hpp>

/*!

*/
//! Array of positions

class CxPositionArray 
{

 public:
  CxPositionArray();
  virtual ~CxPositionArray();
  std::vector<CxPosition> getPositions() const {return (m_positions);};
  bool addPosition(const CxPosition & n_Pos);
	virtual std::vector<CxPosition> getAllAroundPos(const CxPosition& mittelPunkt, float distance )const;
 protected:
  int findPosition(const CxPosition &nPos) const;

 private:
  std::vector<CxPosition> m_positions;



};


#endif
