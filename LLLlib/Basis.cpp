/*! 
  \file Basis.cpp
  
  \brief Implementation of class Basis
  \author Karel Vyborny
*/

/*!
  \class Basis
  \brief A concept for basis, i.e. set of vectors. There are two ways to use it
  
  (1) 'implemented basis' 
  derive a class and specify _what_ are the elements of the basis
  (e.g. yosBasis is a basis consisting of yosStates)

  (2) 'linear combinations' 
  an object of type Basis may be a basis consisting 
  of linear combinations of vectors of some other (possibly implemented)
  basis (like yosBasis).
  \author Karel Vyborny
  \date fall 2002
*/

#include<Basis.hpp>
#include <ERRORS.h>

/*! \fn Basis::Basis(Basis *new_basis,double *new_coef,int new_n)
\brief Construct a basis consisting 
  of linear combinations of vectors of some other basis.

\param int new_n The new basis will have new_n elements, 
\param double *new_coef each element
                 is given by new_basis->dimension() coeficients which 
		 are sequentially stored in new_coef. 
\param new_basis The coeficients
		 are components of the vectors (of the newly 
		 constructed basis) with respect to new_basis.
*/
Basis::Basis(Basis *new_basis,double *new_coef,int new_n):
basis_type (YOS_BASIS),
coef(0)
{
  myBasis=new_basis;
  dim=new_n;
  m_coef.resize(new_n);
 
  if(new_basis != NULL)
    {
    	basis_type = new_basis->getBasisType();
    	
      coef=allocK<double>(myBasis->dimension()*dim);
      for(int i=0;i<dim;i++)
      {
		for(int j=0;j<myBasis->dimension();j++)
		{
	  		coef[i*myBasis->dimension()+j]=new_coef[i*myBasis->dimension()+j];
	  		m_coef.push_back(coef[i*myBasis->dimension()+j]);
      		basis_type=LINEAR_COMBINATION;
		}
      }
    }
}

 
Basis::Basis(const Basis &rhs):
		dim(rhs.dimension()),
		myBasis(0),
basis_type(rhs.getBasisType()),
coef(0),
m_coef()
{
	if (&rhs != this)
	{
		dim = rhs.dimension();	
		coef = rhs.getAllCoefs();
		m_coef = rhs.getAllVCoefs();
		for (int m_index = 0; m_index < dim; ++m_index) {
			
		}
	}
}
Basis::~Basis()
{
  switch(basis_type){
  case LINEAR_COMBINATION:
    delete[] coef;
    break;
  case YOS_BASIS:
    break;
  default:
    throw CxErrors("Basis::~Basis(). Unimplemented type of basis.\n");
  }
}
/*
Basis Basis::operator = (const Basis &rhs)
{
	if (this != &rhs)
	{
		dim = rhs.dimension();
		for(int j=0;j<dim;j++)
		{
	  		coef[i*myBasis->dimension()+j]=new_coef[i*myBasis->dimension()+j];
		}
	}
	return (*this);	
	
}
*/
/*! \fn long int Basis::dimension() const 
  \returns Number of elements of the basis. */
unsigned long int Basis::dimension() const
{
  return dim;
}


void Basis::setDimension(unsigned long  n_dim)
{
	if (n_dim < 0)
	{
		throw CxBadValueError(__FILE__,__LINE__,"Negative dimension in Basis");
	}
	
	dim = n_dim;
	return;
}

/*! \fn double Basis::getCoef(int n_state,int i)
  \brief i-th component of the n_state-th vector of the basis. 
  This makes sense only in Basis and not in derived classes (like yosState).
 */
double Basis::getCoef(int n_state,int i) 
{
  return coef[n_state*myBasis->dimension()+i];
}

/*! \fn double Basis::norm(int i_st)
  Norm of the i_st-th state.
  This makes sense only in Basis and not in derived classes (like yosState).
 */
double Basis::norm(int i_st)
{
  double x=0;
  for(int i=0;i<dim;i++)
    x+=coef[i_st*myBasis->dimension()+i];
  return x;
}

/*! \fn const State *Basis::operator[](long int i)
  Access to the i-th state of the basis. This makes sense only in derived
  classes (overwrite this). */
/*
const State *Basis::operator[](long int i)
{
  throw CxErrors("Basis::operator[]()."
		 "This easy-to-use way is not implemented in Basis."
		 "Do not be lazy and use Basis::getCoef()"); 
	   return 0;
}
*/
/*! \fn BASIS_TYPE Basis::getBasisType()
  \brief Returns the type of basis (see Basis.hpp). */
BASIS_TYPE Basis::getBasisType() const  
{
	return basis_type;
}

/*! \fn Basis *Basis::getMyBasis()
  \brief Returns pointer to basis in which the vectors of myself are given.
  Should be set to NULL in derived classes. */
Basis *Basis::getMyBasis() const
{
	return myBasis;
}

/*! \fn Basis *Basis::getEndBasis() 
  \brief If myBasis is not an implemented basis but still a basis consisting
  of linear combinations, resolve the recursion and find the implemented 
  basis at the end of the chain. */
Basis *Basis::getEndBasis() 
{
  Basis *tmp;
  if(myBasis==NULL) return this;
  tmp=myBasis;
  while(tmp->getMyBasis()!=NULL)
  {
  	 tmp=tmp->getMyBasis();
  }
  return tmp;
}



double * Basis::getAllCoefs() const
{
	return (coef);	
}
//! Accessor to vector member
std::vector <double> Basis::getAllVCoefs() const
{
	return (m_coef);	
}

