#include <utils/CFileParser.h>
#include <utils/logger.hpp>
#include <sstream>

CFileParser::CFileParser()
{
       
}

CFileParser::~CFileParser(void)
{ 

}


CFileParser * CFileParser::m_instance = 0;

CFileParser * CFileParser::getInstance()
{
	if (!m_instance)
	{
		m_instance = new CFileParser();
	}
	return (m_instance);
}
bool CFileParser::hasNextLine(void) const
{
	if (m_fileContent.size()> 0)
	{
		return true;
	}
	else 
	{
		return (false);
	}

}

bool CFileParser::readFile(const std::string& fileName)
{


	m_fileContent.clear();
	glLogger.info("Try to read %s", fileName.c_str());
	std::ifstream inFile;

	try {
	  inFile.open(fileName.c_str(), std::ios_base::in);
	}
	catch (std::exception &e)
	  {
	    glLogger.error("Could not open file");
	    return (false);
	  }
	if (!inFile.is_open())
	{
	  glLogger.warning("Could not open file");
		return false;
	}
	std::string yaString("");
	glLogger.debug("File is open now");
	while (std::getline(inFile,yaString))
	{
		/*
		char  tempChar[255];
		
		inFile.getline(tempChar, 255);


		std::string yaString(tempChar);
		*/

		glLogger.debug("CFileParser read line (%s)", yaString.c_str());
		if (yaString != std::string("") ||!yaString.empty())
		  {
		    m_fileContent.push_back(CLine(yaString));
		  }
	}
	inFile.close();
	return (true);
}

int CFileParser::readFile(const char * fileName)
{
	std::string fileN(fileName);
	return (readFile(fileN));
}

CLine CFileParser::getNextLine(void)
{
	
	CLine retVal = m_fileContent[0];
	m_fileContent.erase(m_fileContent.begin()); 
	return (retVal);
}

