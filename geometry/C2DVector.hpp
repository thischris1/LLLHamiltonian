#ifndef C2DVECTOR_H
#define C2DVECTOR_H
#include <vector>
#include <exception>

using namespace std;

template <class T>
class C2DVector
{
public:
   C2DVector():m_dimRow(0), m_dimCol(0){;}
   C2DVector(int nRow, int nCol) :m_dimRow(nRow), m_dimCol(nCol){

      for (int i=0; i < nRow; i++){
         vector<T> x(nCol);
	 //         int y = x.size();
         m_2DVector.push_back(x);
      }
   }
   void SetAt(int nRow, int nCol, const T& value)  {
     if((unsigned int) nRow >= m_dimRow || (unsigned int)nCol >= m_dimCol)
	//         throw out_of_range("Array out of bound");
	return;
      else
         m_2DVector[nRow][nCol] = value;
   }
  //! Increases value at nRow, nCol by value
  void IncAt(int nRow, int nCol, const T& value) 
  {
    T newVal = GetAt(nRow, nCol) + value;
    SetAt(nRow, nCol, newVal);
    return;
  }
  
   T GetAt( int nRow, int nCol) const {
     if((unsigned int) nRow >= m_dimRow || (unsigned int)nCol >= m_dimCol)
	//         throw out_of_range("Array out of bound");
	return (-42);
      else
         return m_2DVector[nRow][nCol];
   }
   void GrowRow(int newSize) {
      if (newSize <= m_dimRow)
         return;
      m_dimRow = newSize;
      for(int i = 0 ; i < newSize - m_dimCol; i++)   {
         vector<int> x(m_dimRow);
         m_2DVector.push_back(x);
      }
   }
   void GrowCol(int newSize) {
      if(newSize <= m_dimCol)
         return;
      m_dimCol = newSize;
      for (int i=0; i <m_dimRow; i++)
         m_2DVector[i].resize(newSize);
   }
   vector<T>& operator[](int x)    {
      return m_2DVector[x];
   }
  int get_xSize(void) const { return (m_dimRow);}
  int get_ySize(void) const { return (m_dimCol);}
private:
   vector< vector <T> > m_2DVector;
   unsigned int m_dimRow;
   unsigned int m_dimCol;
};

#endif
