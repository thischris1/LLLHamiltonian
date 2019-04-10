#include <geometry/CImpurityArrayValueProvider.h>
CImpurityArrayValueProvider* CImpurityArrayValueProvider::m_instance = 0;


CImpurityArrayValueProvider::~CImpurityArrayValueProvider()
{
}
CImpurityArrayValueProvider::CImpurityArrayValueProvider():
m_array(CImpurityArray())
{
}

bool CImpurityArrayValueProvider::setValue(CImpurityArray &n_array)
{
	m_array = n_array;	
	return (true);
}

double CImpurityArrayValueProvider::getValue(CxPosition &pos) const
{
	
	return (m_array.getPotential(pos));
	
}
CImpurityArrayValueProvider * CImpurityArrayValueProvider::getInstance()
{
	if (m_instance == 0)
	{
		m_instance = new CImpurityArrayValueProvider();
	}
	return m_instance;
}

