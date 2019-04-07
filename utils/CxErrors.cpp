//
/*!
	 \file CxErrors.cpp
 	 \brief The file CxErrors.cpp contains the implementation of the class CxErrors, CxOutOfBoundsError, CxIndextoSmallError, CxIndextoLargeError
         \author Christian Mueller
 	 \date 07 May 03
 */ 

//
//    System Include Files
//

#include <stdio.h>
//
//    Application Include Files
//

#include <utils/CxErrors.hpp>
#include <iostream>
#include <string>
#include <sstream>
//
//    Static Class Variables
//

/*!
	\fn   CxErrors::CxErrors()
	\author Christian Mueller
	\date  07 May 03
 	\brief Standard Constructor
 */


CxErrors::CxErrors()
  : errMsg(std::string(""))
{
    // Test of Parameter


    
    
}

/*!
	\fn   CxErrors::CxErrors(const char * msg)
	\author Christian Mueller
	\date  07 May 03
	\param  msg Message to be displayed
 	\brief 
 */
CxErrors::CxErrors(const char * msg)
{
    if (msg == 0)
    {
        return;
    }
    errMsg = msg;

    
   
}

/*!
  \fn CxErrors::CxErrors(const char * file, int line)
  \param file the file where the error occured (use __FILE__)
  \param line the line where the error occured (use __LINE__)
  \brief builds a short error description (file, line)

 */

CxErrors::CxErrors(const char * file, int line)
{
    if (file != 0)
    {
      std::stringstream mStream;
      mStream << file;
      mStream << " " << line;
      mStream >> errMsg;
    
    }
}



/*!
 \fn CxErrors::CxErrors(const char * message, const char * file,  int line)
 \param message An extra message to be displayed
 \param file The file use __FILE__
 \param line Number where error occured. Use __LINE__ from the calling code normally
 \brief allows to give an error message in the constructor
*/

CxErrors::CxErrors(const char * message, 
		   const char * file, 
		   int line)
{
	std::stringstream mStream;
	if (file != 0)
	{
	
		mStream << "From file " << file << "line No [" << line;
	}
	else 
	{
		mStream << "From unknown file " <<  "line No [" << line;
	}

	if (0 != message) 
	{
		mStream << " [" << message <<"]";
  	}

	mStream >> errMsg;
}

/*!
	\fn void  CxErrors::setMessage(char * newMessage)
	\author Christian Mueller
	\date  21 May 03
	\param  newMessage
 	\brief Sets the error message
 */
/*
	Pre	: 
	Post	: 

 */ 

void CxErrors::setMessage(char * newMessage)
{
    // Test of Parameter
    if (newMessage == 0)
    {
        return;
        
    }

    errMsg = newMessage;


}





/*!
	\fn   CxErrors::~CxErrors()
	\author Christian Mueller
	\date  07 May 03
        
 	\brief The destructor 
 */
/*
	Pre	: 
	Post	: 

 */ 

CxErrors::~CxErrors()
{
    // Do nothing 
  print();


}

CxErrors::CxErrors(const std::string & n_message)
  : errMsg(n_message)
{
}

/*!
	\fn virtual void  CxErrors::print(void)
	\author Christian Mueller
	\date  07 May 03
 	\brief Prints the error message. To be overloaded.
 */


void CxErrors::print(void)
{
    // Test of Parameter


        //! \todo Logging here 
  std::cerr << errMsg.c_str() <<"\n";
        

    
}

const char * CxErrors::getMessage(void) const
{
  return errMsg.c_str();
}
