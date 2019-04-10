#ifndef CIMPURITYARRAYVALUEPROVIDER_H_
#define CIMPURITYARRAYVALUEPROVIDER_H_

#include <utils/IValueProvider.h>
#include <geometry/CImpurityArray.h>
class CImpurityArrayValueProvider : public IValueProvider
{
	
	
public:
	virtual ~CImpurityArrayValueProvider();
	static CImpurityArrayValueProvider * getInstance();
	virtual double getValue(CxPosition &) const;
	bool setValue(CImpurityArray &n_array);
protected:
	CImpurityArrayValueProvider();
	
	
	private:
	CImpurityArray m_array;
	static CImpurityArrayValueProvider * m_instance;
};

#endif /*CIMPURITYARRAYVALUEPROVIDER_H_*/
