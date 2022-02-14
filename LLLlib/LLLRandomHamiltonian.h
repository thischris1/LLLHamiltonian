#ifndef LLLRANDOMHAMILTONIAN_H_
#define LLLRANDOMHAMILTONIAN_H_

#include "LLLhamiltonian.hpp"

class LLLRandomHamiltonian : public LLLhamiltonian
{
public:
	LLLRandomHamiltonian(yosBasis *new_basis,int new_type,double bli_new,
		 double a_new, double b_new);
	virtual ~LLLRandomHamiltonian();
};

#endif /*LLLRANDOMHAMILTONIAN_H_*/
