#include <utils/CxIndextoSmallError.hpp>
#include <string>
#include <sstream>


CxIndextoSmallError::CxIndextoSmallError(const int Index)
  :CxOutOfBoundsError(createErrorMessage(Index))
{
   
    
}
/*!
	\fn   CxIndextoSmallError::CxIndextoSmallError()
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

CxIndextoSmallError::CxIndextoSmallError(const char * file , const int lineNo) : CxOutOfBoundsError(file, lineNo)
{
    

}
/*!
\fn CxIndextoSmallError::~CxIndextoSmallError()
\brief destructor

*/

CxIndextoSmallError::~CxIndextoSmallError()
{
  print();
}



std::string  CxIndextoSmallError::createErrorMessage( const int index)
{
  std::stringstream tempStream;
  tempStream << "Index ("<<index<<") is too small";
  std::string retVal;
  tempStream >> retVal;
  return retVal;
}
