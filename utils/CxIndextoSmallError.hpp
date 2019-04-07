#ifndef CxIndextoSmallError_hpp
#define CxIndextoSmallError_hpp

#include <utils/CxErrors.hpp>
#include <string>
#include <utils/CxOutOfBoundsError.hpp>

/*!
  \class CxIndextoSmallError
  \brief thrown for negative aray indices
*/

class CxIndextoSmallError : public CxOutOfBoundsError
{

public :
  CxIndextoSmallError(const int Index);
  CxIndextoSmallError(const char * file, const int lineNo);
  virtual ~CxIndextoSmallError();
  static std::string  createErrorMessage( const int index);

    
private:
    int itsIndex;
    
    char* errMsg;
};










#endif
