#ifndef CLINE_H
#define CLINE_H
#include <string>
#include <vector>


class CLine
{
public:
	CLine(void);
	virtual ~CLine(void);
private:
	// The line read from the parser
	std::string m_line;
public:
	CLine(std::string newLine);
	bool setLine(std::string n_line);
	bool isCommentLine(void);
	std::vector<int> getIntegers(void) const;
	// returns the whole line
	std::string getLine(bool withComment = false) const;
private:
	char m_commentChar;
protected:
	std::string stripComment(void) const;
public:
	std::vector<float> getFloats(void) const;
	std::string getString(void) const;
};
#endif
