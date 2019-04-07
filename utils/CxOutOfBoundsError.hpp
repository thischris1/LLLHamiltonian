#ifndef CXOUTOFBOUNDSERROR_HPP
#define CXOUTOFBOUNDSERROR_HPP
#include <utils/CxErrors.hpp>



/*!
  \class CxOutOfBoundsError : public CxErrors
  \brief Class to be thrown for an out of bounds access



 */
class CxOutOfBoundsError : public CxErrors
{
public:
    
  CxOutOfBoundsError();
  CxOutOfBoundsError( char const * file, int lineNo): CxErrors(file, lineNo){};
  CxOutOfBoundsError( const std::string & n_message) : CxErrors(n_message){};
  virtual ~CxOutOfBoundsError();


  
};







#endif
