#include <utils/CxIndextoLargeError.hpp>
/*!
  \fnCxIndextoLargeError::~CxIndextoLargeError(
  \brief Destructor

*/
CxIndextoLargeError::~CxIndextoLargeError()
{

}
CxIndextoLargeError::CxIndextoLargeError(const int Index)
        :CxOutOfBoundsError("",0)
{
 
    
}
/*!
	\fn   CxIndextoLargeError::CxIndextoLargeError()
	\author Christian Mueller
	\date  14 May 03
	\return 
	\param 
 	\brief 
 */
/*
	Pre	: 
	Post	: 

 */ 

CxIndextoLargeError::CxIndextoLargeError(const char * file , const  int lineNo) : CxOutOfBoundsError(file, lineNo)
{
    

}




