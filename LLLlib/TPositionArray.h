#ifndef TPOSITIONARRAY_H
#define  TPOSITIONARRAY_H
#include <vector>
template <class T> 
class TPositionArray
{
 public:
  TPositionArray() {};
  
  virtual ~TPositionArray() 
    {
      m_positions.clear();
    };
  
  std::vector<T> getPositions() const {return (m_positions);};


  //! Add a position

  bool addPosition(const  T  & n_Pos)
  {
      if (findPosition(n_Pos) !=-1)
    {
      return (false);
    }
  m_positions.push_back(n_Pos);
  return (true);
  }


  unsigned int size() const { return m_positions.size();};


  void clear(void) { m_positions.clear();};


 protected:
  
  int findPosition(const T  &nPos) const
  {
    int retVal = -1;
    
    for (unsigned int index = 0; index < m_positions.size(); index++)
      
      {
	T  m_pos = m_positions[index];
	if (nPos.compare(m_pos))
	{
	  retVal = index;
	  break;
	}

    }
  return (retVal);
}





 private:
  std::vector<T> m_positions;



};

#endif
