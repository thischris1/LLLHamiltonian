//
/*!
	 \file logger.hpp
 
 	 \brief The file "logger.hpp" contains the declaration of the class Logger 

  	  \author Christian Mueller
          
          \date 07 Feb 03
          
*/ 

//
//    System Include Files
//

#ifndef LOGGER_HPP
#define LOGGER_HPP


#include <fstream>
#include <stdio.h>
#include <cstdarg>  //for variable argument lists

//
//    Application Include Files
//
#include <LLLlib/yosState.hpp>
#include <geometry/CxImpurity.hpp>


//
//    Global Variables
//
#define ERROR(msg) glLogger.error(msg);

#define WARNING(msg) glLogger.warning(msg);

#define INFO(msg) glLogger.info(msg);

#define DEBUG(msg)  glLogger.debug(msg);
#define ERROR_PARAM(msg,param) glLogger.error(msg,param)

/*!
	\class Logger
	\brief The class Logger 
	\author Christian Mueller
	\date 07 Feb 03
        This class is intended to be used with all other classes as a common log-output interface. It has 4 Loglevels (for now) INFO, DEBUG, ERROR, WARNING

*/
using namespace std;

typedef enum
{
    ERROR,
    WARNING,
    INFO,
    DEBUG,
    MAXLOGLEVEL
} Loglevel; //! Determines the loglevel




    

class Logger 
{
public:
/*! 
  \fn  Logger()
  \brief Constructor
        
*/
    Logger();
    
    
    /*! 
 	\fn     int init( const char *,  constc char *)
	\return    0 if successfully initialized the Logger (i.e. all files are there and are writeable.
        
	\brief Initializes the Logger
        
    */
    /*! 
 	\fn      Logger(const char* fileName)
	\return    
	\param const char* fileName    
	\brief 

 */ 
    explicit Logger(const char* fileName, bool fileOnlyLog = false);



  //!  Set the Logfile after construction

  int setFileName(const char *);
  
  int  init( const char *, const char*);
    
    /*! 
      \fn     void warning(const char*)
      \return    nothing
      \brief Writes a warning to the screen and/or logfile
      
    */ 


  //! sets the Log-level.

  void setLogLevel(const int newLogLevel);
  //! Set loglevel with enum value
  void setLogLevel(const Loglevel);

  //! Set loglevel with char
  void setLogLevel(const char * c_logLevel);

//! Set loglevel with string
  void setLogLevel(const std::string &logLevel);
  //! Disables console logging
  void disableConsole(void);
  
  //! Enables console logging
  void enableConsole(void);
/*! 
 	\fn    LogLevel getLogLevel(void) const
	\return The Loglevel      
	\brief Returns the Loglevel of the class

 */ 
  Loglevel getLogLevel(void) const;
 

/*! 
 	\fn     void warning(const char*)
	\return    
	\param const char* text containing format specifiers like in printf   
	\param ... arbitrary number of additional arguments to be output in format
	given in text
	\brief 

 */ 
  void  warning(const char* ...);

  /*!
    \fn void warning(const yosState & )
  */
  void warning(const yosState & );
  /*!
    \fn void warning(const CxImpurity & )
  */
  void warning(const CxImpurity & );

/*! 
 	\fn     void error(const char*)
	\return    
	\param const char* text containing format specifiers like in printf   
	\param ... arbitrary number of additional arguments to be output in format
	given in text
	\brief 

 */ 
    void  error(const char* ...);

  /*!
    \fn void error(const yosState)
  */

  void error(const yosState &);

  void error (const std::string &);

/*!
    \fn void error(const CxImpurity &)
  */

  void error(const CxImpurity &);

/*! 
 	\fn     void info(const char * ...)
        \param ... arbitrary number of additional arguments to be output in format given in text
	\brief Logs an arbitary list of values (printf like)

 */ 
    void  info(const char * ...);


  /*!
    \fn void info(const yosState)

  */
  
  void info(const yosState &);

 /*!
    \fn void info(const CxImpurity &)

  */
  
  void info(const CxImpurity &);

/*! 
 	\fn     void debug(const char *, ...)
	\return    
	\param const char* text containing format specifiers like in printf   
	\param ... arbitrary number of additional arguments to be output in format
	given in text
	\brief

 */ 
    void  debug(const char * ...);

/*!
    \fn void debug(const yosState &)

  */
  void debug(const yosState &);


  /*!
    \fn void debug(const CxImpurity &)

  */
  void debug(const CxImpurity &);  
  /*!
    \fn     : ~Logger
    \brief Destructor
  */
        
    ~Logger();
  //! Switches the logging header on  
  void enableHeader(void);
  //! Switches the logging header on  
  void disableHeader(void);
protected:

  static const unsigned int MESSAGE_LEN_GUESS = 800;
  static const unsigned int MESSAGE_MAX_LEN_WITHOUT_PARS = 200;

    /*! 
 	\fn     void showLogMessage(const int &, const char *, ...)
	\return    
	\param const char * format string containing text and printf-like format specifiers
	\param ... variable number of arguments (like in printf)
	\brief 
    */ 
  void  showLogMessage(const int &, const char *, va_list);

  /*!
    \fn  void showLogMessage(const yosState &, const LogLevel theLevel) const
  */
  void showLogMessage(const yosState &, 
		      const Loglevel & theLevel) const;


  /*!
    \fn  void showLogMessage(const CxImpurity &, const LogLevel theLevel) const
  */
  void showLogMessage(const CxImpurity &, 
		      const Loglevel & theLevel) const;

/*!
  \fn   void printLogMessage(Loglevel logLevel) const;
  \brief small helper function, often Overloaded
 */
  
    void printLogMessage(Loglevel logLevel) const;

    void printLogMessage(const char* message) const;

    void printLogMessage(const int value) const;
    
  private:
    /*!
        
      \fn     : Logger(const Logger& rOrig);
      \brief Copy Constructor
      
      \param rOrig the value to be copied. Nothing happens here, it is not intended to be copied!;
    */
    
    
    Logger(const Logger& rOrig);
    
/*!
  \fn     : operator = 
  \brief Assignment Operator
*/
    Logger& operator = (const Logger& rOrig);

/*! 
 	\fn     char * getTimeString(int &, char * )
	\return The current time as a string in HH:MM:SS.MS format.
        The length of the string is returned in the integer parameter
	
 */ 
   std::string getTimeString(void ) const;





    
//! Member variables 
    
  Loglevel myLogLevel;
    

  ofstream *logFile;
  bool hasHeader;
  bool fileOnly;
  };


// one global instance of Logger class declared in logger.cpp
// may be used in any cpp-file that includes <logger.hpp> 
extern Logger glLogger;



#endif //LOGGER_HPP
