#ifndef CONSTANTVALUEPROVIDER_H_
#define CONSTANTVALUEPROVIDER_H_

#include "IValueProvider.h"

class ConstantValueProvider : public IValueProvider
{
public:
	
	virtual ~ConstantValueProvider();
	static ConstantValueProvider * getInstance(){return (m_instance);};
	virtual double getValue(CxPosition &) const;
	bool setValue(double n_value){m_value=n_value; return (true);};

	private:
	ConstantValueProvider();

	static ConstantValueProvider * m_instance;
	double m_value;
};

#endif /*CONSTANTVALUEPROVIDER_H_*/
