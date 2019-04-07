#ifndef CXERRORS_HH
#define CXERRORS_HH



/*!
	\file CxErrors.hpp
	\brief The file  contains the declaration of the class CxErrors,   CxOutOfBoundsError, CxIndextoSmallError, CxIndextoLargeError, CxFileNotFoundError
	\author Christian Mueller
	\date 07 May 03

	
*/ 
 
//
//    System Include Files
//
#include <string>
/*!
  \class  CxErrors
  \brief Abstract basis class for error handling. It carries a message that is printed at request
*/

class CxErrors
{
public:
  CxErrors();
  CxErrors(const char * file,  int lineNo);
  CxErrors(const char * msg);
  CxErrors(const char * message, const char * file,  int lineNo);
  CxErrors(const std::string & n_message); 
  void setMessage(char *newMsg);
  const char * getMessage(void) const;

  void print(void);
  virtual ~CxErrors();   
protected:
 

private:
    
  std::string errMsg;
    

};

#endif // CXERRORS_HH
