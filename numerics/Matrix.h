#ifndef MATRIX_H_
#define MATRIX_H_
#include <LLLlib/LLLlib.h>
#include <LLLlib/Types.hpp>
class Matrix
{
public:
	Matrix(int rows, int col);
	virtual ~Matrix();
	int getRow(void)const {return (m_rows);};
	int getCols(void) const {return (m_cols);};
	virtual CMDCOMPLEX getElement(const int row, const int col)=0;
	virtual bool setElement(const CMDCOMPLEX  n_value, const int row, const int col)=0;
private:
	int m_rows;
	int m_cols;
	
	// Methods wo implementation
	Matrix();

};

#endif /*MATRIX_H_*/
