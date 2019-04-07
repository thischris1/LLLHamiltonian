#include "ConstantValueProvider.h"

ConstantValueProvider * ConstantValueProvider::m_instance = new ConstantValueProvider();


ConstantValueProvider::ConstantValueProvider():
m_value(1.0)
{
}

ConstantValueProvider::~ConstantValueProvider()
{
}



double ConstantValueProvider::getValue(CxPosition &) const
{
	return (m_value);
}

/*
ConstantValueProvider* ConstantValueProvider:getInstance()
{
	if (!m_instance)
	{
		m_instance = new ConstantValueProvider();
	}
	return (m_instance);
}
*/

