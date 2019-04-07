#ifndef CXBADVALUEERROR_HPP
#define CXBADVALUEERROR_HPP
#include "CxErrors.hpp"



/*!
  \class CxBadValueError
*/

class CxBadValueError:public CxErrors
{
public:
  
  CxBadValueError();
  CxBadValueError(char const * file, int lineNo);
  CxBadValueError(char const * file, int lineNo, const char *message);
  virtual ~CxBadValueError();
  void handle();



};




#endif
