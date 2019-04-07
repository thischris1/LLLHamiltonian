#ifndef COMPLEXSQUAREMATRIX_H_
#define COMPLEXSQUAREMATRIX_H_
#include <numerics/SquareMatrix.h>
#include <LLLlib/Types.hpp>
#include <vector>
class ComplexSquareMatrix :public SquareMatrix
{
public:
	ComplexSquareMatrix(int rank);
	virtual ~ComplexSquareMatrix();
	
private:
	ComplexSquareMatrix();
	// members
	std::vector<CMDCOMPLEX> m_values;
	std::vector<int> m_rowIndex;
	std::vector<int> m_colIndex;
};

#endif /*COMPLEXSQUAREMATRIX_H_*/
