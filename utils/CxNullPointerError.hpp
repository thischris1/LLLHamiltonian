#ifndef CxNullPointerError_hpp
#define CxNullPointerError_hpp
#include <utils/CxErrors.hpp>
/*!

  \class CxNullPointerError : public CxErrors
  \brief to be thrown when a null pointer was passed in a method
  
 */

class CxNullPointerError : public CxErrors
{
public:
  CxNullPointerError():CxErrors("NullPointerError"){};
  CxNullPointerError(const char * file, const int lineNo): CxErrors("Nullpointererror",file, lineNo){};
    CxNullPointerError(const char * message,const char * file,int  lineNo) 
      : CxErrors(message, file, lineNo){};
     ~CxNullPointerError(void);



};


#endif
