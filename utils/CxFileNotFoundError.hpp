#ifndef CxFileNotFoundError_hpp
#define CxFileNotFoundError_hpp
#include <utils/CxErrors.hpp>
/*!

  \class CxFileNotFoundError : public CxErrors
  \brief to be thrown when a null pointer was passed in a method
  
 */

class CxFileNotFoundError : public CxErrors
{
public:
    CxFileNotFoundError();
  CxFileNotFoundError(const char * file, int lineNo): CxErrors("FILE not found", file, lineNo){};
    CxFileNotFoundError(const char * message,
        const char * file,
			int lineNo): CxErrors(message, file, lineNo){};
    virtual ~CxFileNotFoundError(void);

};

#endif
