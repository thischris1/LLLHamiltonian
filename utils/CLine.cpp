#include <utils/CLine.h>
#include <sstream>
#include <utils/logger.hpp>


CLine::CLine(void)
: m_commentChar('#')
{
	m_line.clear();
}

CLine::~CLine(void)
{
	m_line.clear();
}

CLine::CLine(std::string newLine):
		 m_line("#"),
		 m_commentChar('#')
{
	setLine(newLine);

}

bool CLine::setLine(std::string n_line)
{
	m_line = n_line;
	return (true);
}

bool CLine::isCommentLine(void)
{
	return false;
}


/*!
	returns the numbers in line as integers
*/
std::vector<int> CLine::getIntegers(void) const
{
	glLogger.debug("Parsing line");
	glLogger.debug(getLine().c_str());
	// Strip spurious comments
	std::string strippedString  = stripComment();
	std::stringstream myStream;

	myStream<< strippedString;
		
	std::vector<int> integerList;	

	int tempVal;
	while (myStream >> tempVal)
	{
		if (!myStream.fail())
		{
			integerList.push_back(tempVal);
		}
	}
	if (!myStream.eof())
	{
		integerList.pop_back();
	}
	return (integerList);
}

// returns the whole line
std::string CLine::getLine(bool withComment) const
{
	if (!withComment)
	{
		return (m_line);
	}
	return (stripComment());
	
}

std::string CLine::stripComment(void) const
{
	return (m_line.substr(0, m_line.find(m_commentChar)));
}

std::vector<float> CLine::getFloats(void) const
{
	std::string strippedString  = stripComment();
	//std::istringstream
glLogger.debug("Parsing line for floats");
glLogger.debug(getLine().c_str());
	std::stringstream myStream;
	std::vector<float> retVal;
	
	myStream.setf(std::ios::dec, std::ios::scientific);
	myStream<< strippedString;
	float tempVal;
	while ( myStream >> tempVal)
	{
		retVal.push_back(tempVal);
	}
	if (!myStream.eof())
	{
		retVal.pop_back();
	}
	return (retVal);
}

std::string CLine::getString(void) const
{
	
	return (stripComment());
}
