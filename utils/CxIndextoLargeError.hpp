#ifndef CxIndextoLargeError_hpp
#define CxIndextoLargeError_hpp
#include <utils/CxOutOfBoundsError.hpp>
#include <string>

/*!
  \class CxIndextoLargeError : public CxOutOfBoundsError
  \brief to be thrown for exceeding array sizes or index reading
*/
class CxIndextoLargeError : public CxOutOfBoundsError
{

public :
  CxIndextoLargeError(const int Index);
  CxIndextoLargeError(const char * file, const int lineNo);
  virtual ~CxIndextoLargeError();
  
    

    
};
    
#endif
