#ifndef CFILEPARSER_H
#define CFILEPARSER_H


#define CFILEPARSER_H
#include <fstream>
#include <vector>
#include <string>
#include <utils/CLine.h>

class CFileParser 
{
public:
	virtual ~CFileParser(void);
public:
	static CFileParser * getInstance(void);
	

private: 
	CFileParser();
	static CFileParser * m_instance;
	

public:
	bool hasNextLine(void) const;
private:
//	std::vector<std::string> m_fileContent;
	std::vector<CLine> m_fileContent;
public:
	bool readFile(const std::string& fileName);
	int readFile(const char * fileName);
private:
public:
	CLine getNextLine(void);
};
#endif
