#include <utils/CxBadValueError.hpp>

/*!

*/

CxBadValueError::CxBadValueError():CxErrors("Bad value error occured")
{

}
/*!

*/
CxBadValueError::CxBadValueError(const char * fileName,
				 int lineNo)
  :CxErrors (fileName, lineNo)
{
    
    
}

CxBadValueError::CxBadValueError(const char * fileName,
				 int lineNo, 
				 const char * message)
  : CxErrors(message, fileName, lineNo)
{

    
}

CxBadValueError::~CxBadValueError()
{


}
