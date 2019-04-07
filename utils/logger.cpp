/*!

	 \file logger.cpp
 	 \brief The file  contains the implementation of the class 
         \author Christian Mueller
 	 \date 07 Feb 03

 */ 

//
//    System Include Files
//
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <malloc.h>
#include <string>
//
//    Application Include Files
//

#include <utils/logger.hpp>


/*

  MACROS

 */
using namespace std;

/*!
	\fn (const char* message)  Logger::warning(const char* text)
	\author Christian Mueller
	\date  07 Feb 03
	\return 
	\param const char* text
 	\brief 

	\todo Check for too long messages in showLogMessage 
	(use something like vasprintf() which would be standard). 
	vasprintf() is not known at alphas.
 */
/*
	Pre	: 
	Post	: 

 */ 

Logger::Logger():
		myLogLevel(DEBUG),
		logFile(0),
		hasHeader(true),
		fileOnly(false)
{

}


Logger::~Logger()
{
	if (logFile != NULL)
	{
		logFile->close();
		delete logFile;
		logFile = 0;
	}
	delete logFile;
	logFile = NULL;
}

/*!
	\fn   Logger::Logger(const char* fileName)
	\author Christian Mueller
	\date  07 Feb 03
	\param  fileName Name of file to write to
	\param  fileOnlyLog set it to true if you dont want any console-logging
 	\brief constructor with a file to write logging to.
 */
/*
	Pre	: 
	Post	: 
Loglevel myLogLevel;


  ofstream *logFile;
  bool hasHeader;
  bool fileOnly;
 */ 

Logger::Logger(const char* fileName,
		bool fileOnlyLog):
				myLogLevel(DEBUG),
				logFile(0),
				hasHeader(true),
				fileOnly(fileOnlyLog)
{
	// Test of Parameter

	setFileName(fileName);


}

/*!
	\fn   int Logger::setFileName(const char* fileName)
	\author Christian Mueller
	\date  07 Feb 03
	\param  fileName The name of the logfile
 	\brief sets the File to which the logger is supposed to write.
	\return 0 if successfully set 1 otherwise
 */

int Logger::setFileName(const char* fileName)
{
	if (fileName != 0 )
	{
		// Make that more comfortable later (i.e. Check for directory, readability etc.
		if (logFile == 0)
		{
			logFile = new ofstream(fileName);
			return 0;
		}
		else
		{
			error("you can not change the name of the file.");
			return 1;
		}
	}
	else
	{
		error("You passed an empty string as logfile."
				"Output is written to screen only");
		return 1;
	}


}


/*!
  \fn void Logger:setLogLevel(const LogLevel newLogLevel)
  \brief sets the Loglevel to newLogLevel
  \param  newLogLevel The new Loglevel

 */

void Logger::setLogLevel(const Loglevel newLogLevel)
{
	cerr << "LOGLEVEL SET TO  newLogLevel "<<  newLogLevel <<"\n";
	error("LOGLEVEL SET (%d)", newLogLevel);
	myLogLevel = newLogLevel;
}




/*!
  \fn void Logger:setLogLevel(const int  newLogLevel)
  \brief sets the Loglevel to newLogLevel
  \param  newLogLevel The new Loglevel

 */

void Logger::setLogLevel(const int newLogLevel)
{
	cerr << "LOGLEVEL SET TO  newLogLevel "<<  newLogLevel <<"\n";
	error("LOGLEVEL SET (%d)", newLogLevel);
	switch (newLogLevel) {
	case (0):
    		myLogLevel= ERROR;
	break;
	case (1):
    		myLogLevel= WARNING;
	break;
	case(2):
    		myLogLevel=INFO;
	break;
	default:
		myLogLevel=INFO;
		break;
	}
}

void Logger::setLogLevel(const std::string &logLevel)
{


	if (std::string("0") == logLevel)
	{
		myLogLevel= ERROR;
		return;
	}
	else if (std::string("1") == logLevel)
	{
		myLogLevel= WARNING;
		return;
	}
	else if (std::string("2") == logLevel)
	{
		myLogLevel= INFO;
		return;
	}
	else if (std::string("3") == logLevel)
	{
		myLogLevel= DEBUG;
		return;
	}
	else
	{
		myLogLevel= ERROR;
		return;
	}

}

void Logger::setLogLevel( const char * c_logLevel)
{
	setLogLevel(std::string(c_logLevel));

}
/*!
  \fn void Logger::disableConsole(void)
 */

void Logger::disableConsole(void)
{
	fileOnly = true;
}

/*!
  \fn void Logger::enableConsole(void)
 */

void Logger::enableConsole(void)
{
	fileOnly = true;
}



/*!
  \fn LogLevel Logger::getLogLevel(void) const
 */



Loglevel Logger::getLogLevel(void) const
{
	return myLogLevel;
}


/*!
	\fn void  Logger::error(const char * text ...)
	\author Christian Mueller
	\date  07 Feb 03
	\return nothing 
	\param  text The Logstring containing formt specifiers like in printf and additional parameters
 	\brief Writes the string to the log-stream (either screen or file or both) if loglevel is ERROR or higher (i.e. always)

 */

void Logger::error(const char * text ...)
{

	if (myLogLevel >= ERROR)
	{
		va_list va;
		va_start (va, text);  // convert ellipsis to variable argument list
		showLogMessage(ERROR, text, va);
		va_end (va);
	}
	else
	{
		return;

	}
	return;

}

/*!
	\fn void  Logger::error(const yosState logState)
	\author Christian Mueller
	\date  07 Feb 03
	\return nothing 
	\param  logState. A yosState to be logged. s.a. error()
 	\brief Writes the string to the log-stream (either screen or file or both) if loglevel is ERROR or higher (i.e. always)

 */
void Logger::error(const yosState & logState)
{
	if (myLogLevel >= ERROR)
	{
		showLogMessage(logState, ERROR);

	}


}

void Logger::error(const std::string & theString)
{
	error(theString.c_str());
}

/*!
	\fn void  Logger::info(const char * text, ...)
	\author Christian Mueller
	\date  07 Feb 03
	\return nothing 
	\param const char * text
 	\brief 
 */

void Logger::info(const char * text ...)
{

	if (myLogLevel >= INFO)
	{
		va_list va;
		va_start (va, text);
		showLogMessage(INFO, text, va);
		va_end (va);
	}

}

/*!
	\fn void  Logger::info(const yosState logState)
	\author Christian Mueller
	\date  07 Feb 03
	\return nothing 
	\param  logState. A yosState to be logged. s.a. error()
 	\brief Writes the string to the log-stream (either screen or file or both) if loglevel is ERROR or higher (i.e. always)

 */
void Logger::info(const yosState & logState)
{
	if (myLogLevel >= INFO)
	{
		showLogMessage(logState, INFO);

	}

}


/*!
	\fn void  Logger::debug(const char * text ...)
	\author Christian Mueller
	\date  07 Feb 03
	\return nothing 
	\param  text ... to be printed
 	\brief 
 */


void Logger::debug(const char * text ...)
{

	// Test of Parameter

	if (myLogLevel >= DEBUG)
	{
		va_list va;
		va_start (va, text);

		showLogMessage(DEBUG, text, va);

		va_end (va);
	}
}

/*!
  \fn void  Logger::debug(const yosState logState &)
  \author Christian Mueller
  \date  07 Feb 03
  \return nothing 
  \param  logState. A yosState to be logged. s.a. error()
  \brief Writes the string to the log-stream (either screen or file or both) if loglevel is ERROR or higher (i.e. always)

 */
void Logger::debug(const yosState & logState)
{
	if (myLogLevel >= DEBUG)
	{
		showLogMessage(logState, DEBUG);

	}

}


/*!
  \fn void Logger::warning(const char* text ...)
 */

void Logger::warning(const char* text ...)
{
	// Test of Parameter
	if (myLogLevel >= WARNING)
	{
		va_list va;
		va_start (va, text);

		//Add a timestamp
		showLogMessage(WARNING, text, va);

		va_end (va);
	}
	return;
}

/*!
  \fn void  Logger::warning(const yosState logState &)
  \author Christian Mueller
  \date  07 Feb 03
  \return nothing 
  \param  logState. A yosState to be logged. s.a. error()
  \brief Writes the string to the log-stream (either screen or file or both) if loglevel is ERROR or higher (i.e. always)

 */
void Logger::warning(const yosState & logState)
{
	if (myLogLevel >= WARNING)
	{
		showLogMessage(logState, WARNING);

	}

}
/*!
  \fn void  Logger::warning(const yosState logState &)
  \author Christian Mueller
  \date  07 Feb 03
  \return nothing 
  \param  logState. A yosState to be logged. s.a. error()
  \brief Writes the string to the log-stream (either screen or file or both) if loglevel is ERROR or higher (i.e. always)

 */
void Logger::warning(const CxImpurity & logState)
{
	if (myLogLevel >= WARNING)
	{
		showLogMessage(logState, WARNING);

	}

}

void Logger::error(const CxImpurity & logState)
{

	showLogMessage(logState, ERROR);


}

void Logger::info(const CxImpurity & logState)
{
	if (myLogLevel >= INFO)
	{
		showLogMessage(logState, INFO);

	}

}

void Logger::debug(const CxImpurity & logState)
{
	if (myLogLevel >= INFO)
	{
		showLogMessage(logState, DEBUG);

	}

}

/*!
	\fn void  Logger::showLogMessage(const int& logLevel, const char * text, va_list va)
	\author Christian Mueller
	\date  07 Feb 03
	\return nothing
        \param logLevel the Loglevel
	\param text containing format specifiers
	\param va variable argument list of additional passed arguments
 	\brief Composes the logmessage and prints it to STDERR and/or file
 */
/*
	Pre	: 
	Post	: 

 */ 
void Logger::showLogMessage(const int & logLevel, const char * text, va_list va)
{

	// Test of Parameter

	if (text == 0)
	{
		return;
	}

	//  complete message out of variable argument list
	// use automatically allocating printf, remember to free in the end !
	char *new_text =new char[MESSAGE_LEN_GUESS];

	/*if( strlen(text) > MESSAGE_MAX_LEN_WITHOUT_PARS)
      {
	cerr << 
	  "Dangerously long string passed to logger (l = " 
	     << strlen(text) << ")\n";
      }
	 */
	vsprintf (new_text, text, va);
	//#endif


	if (hasHeader)
	{
		const char * timeBuffer = getTimeString().c_str();

		// get and Convert time
		printLogMessage(timeBuffer);
	}

	printLogMessage((Loglevel)logLevel);
	printLogMessage(new_text);
	printLogMessage("\n");

	delete[] new_text;
}


/*!
  \fn void showLogMessage (const yosState & logState) const

 */
void Logger::showLogMessage (const yosState & logState, 
		const Loglevel & logLevel)const
{
	if (logState.getNe() <= 0)
	{
		return;
	}



	std::string timeBuffer = getTimeString().c_str();



	// get and Convert time

	if (hasHeader)
	{

		printLogMessage(timeBuffer.c_str());
		printLogMessage(logLevel);
	}

	printLogMessage( " | ");
	bool spin = logState.hasSpin();
	for (int iIdx = 0;
			iIdx < logState.getNe();
			iIdx++)
	{
		printLogMessage(logState.j(iIdx));
		if (spin)
		{
			if (logState.spin(iIdx) == 0)
			{
				printLogMessage("-") ;
			}
			else
			{
				printLogMessage("+") ;
			}
		}

	}

	printLogMessage(" > \n");

}


/*!
  \fn void showLogMessage (const CxImpurity & logState) const

 */
void Logger::showLogMessage (const CxImpurity & logState, 
		const Loglevel & logLevel)const
{

	if (hasHeader)
	{
		const char * timeBuffer = getTimeString().c_str();



		// get and Convert time

		printLogMessage(timeBuffer);
		delete [] timeBuffer;
		printLogMessage(logLevel);
	}

	char * logMessage = new char[100];
	sprintf (logMessage, "Impurity at x = (%f), y=(%f), sigmaX=(%f), sigmaY(%f)"
			" Strength = (%f) \n",
			logState.getXPosition(),
			logState.getYPosition(),
			logState.getSigmaX(),
			logState.getSigmaY(),
			logState.getStrength());

	printLogMessage( logMessage);
	delete [] logMessage;



}







/*!
	\fn void  Logger::printLogMessage(Loglevel logLevel) const
	\author Christian Mueller
	\date  15 Apr 03
	\param  logLevel 
 	\brief Prints the Loglevel as string
 */
/*
	Pre	: 
	Post	: 

 */ 

void Logger::printLogMessage(Loglevel logLevel) const
{
	if (hasHeader)
	{
		switch
		(logLevel)
		{
		case INFO:
		{
			if (!fileOnly)
			{

				cerr << "INFO: ";
			}

			if (logFile)
			{
				(*logFile) << "INFO: ";
			}
			break;
		}

		case ERROR:
		{
			if (!fileOnly)
			{
				cerr << "ERROR: ";
			}

			if (logFile)
			{
				(*logFile) << "ERROR: ";
			}
			break;
		}

		case WARNING:
		{
			if (!fileOnly)
			{
				cerr << "WARNING: ";
			}
			if (logFile)
			{
				(*logFile) << "WARNING: ";
			}
			break;
		}
		default:
		{
			if (!fileOnly)
			{
				cerr << "DEBUG: ";
			}
			if (logFile)
			{
				(*logFile) << "DEBUG: ";
			}


			break;
		}
		}
	}
}
/*!
	\fn void  Logger::printLogMessage(const int value) const
	\author Christian Mueller
	\date  15 Apr 03
	\param  value The index to be printed 
 	\brief Prints an integer, small helper function to be used with other logging functions
 */


void Logger::printLogMessage(const int value) const
{

	if (logFile)
	{
		(*logFile) << value << " ";
	}
	if (!fileOnly)
	{
		cerr << value << " ";
	}

}





/*!
	\fn void  Logger::printLogMessage(const char* message) const
	\author Christian Mueller
	\date  15 Apr 03
	\param  message The message to be printed
 	\brief writes message to screen and file

        Small helper function
 */


void Logger::printLogMessage(const char* message) const
{
	if (message == 0)
	{
		return;

	}

	if (logFile)
	{
		(*logFile) << message << " ";
	}
	if (!fileOnly)
	{
		cerr <<message << " ";
	}



}





/*!
	\fn void  Logger::getTimeString(int & stringLength, char * retVal) const 
	\author Christian Mueller
	\date  07 Feb 03
        \param  retVal The string with the date
	\param stringLength returns the length of the string returned
 	\brief Reads the time from the system and converts it into a string. Non portable to NON_UNIX machine (I assume)
 */


std::string Logger::getTimeString(void) const 
{
	// Test of Parameter
	/*
    char *retVal = new  char[20];

    time_t t = time(0);
    strftime (retVal, 20, "%X", localtime (&t));
	 */
	//    std::string timeString(retVal);
	std::string timeString("");
	//    delete [] retVal;

	return timeString;


}

void setLogLevel( const int newLoglevel)
{

	switch (newLoglevel)
	{

	case(0):
    	  {
		setLogLevel(ERROR);
		break;
    	  }
	case (1):
    	  {
		setLogLevel(WARNING);
		break;
    	  }
	case (2):
    	  {
		setLogLevel(INFO);
		break;
    	  }
	case(3):
    	  {
		setLogLevel(DEBUG);
		break;
    	  }
	default:
	{
		//warning("There is no such loglevel");
		setLogLevel(DEBUG);
	}
	return;

	}



}

void Logger::disableHeader(void)
{
	hasHeader = false;
}

void Logger::enableHeader(void)
{
	hasHeader = true;
}


// one global accessible instance of Logger
// made accessible through external declaration in logger.hpp
Logger glLogger;
