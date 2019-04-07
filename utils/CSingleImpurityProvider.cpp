#include <utils/CSingleImpurityProvider.h>

CSingleImpurityProvider* CSingleImpurityProvider::m_instance = new CSingleImpurityProvider();

CSingleImpurityProvider::~CSingleImpurityProvider()
{
}
CSingleImpurityProvider::CSingleImpurityProvider():
m_imp(CxImpurity(0.5,0.5,0.1,0.1,-0.3))
{
}

bool CSingleImpurityProvider::setValue(CxImpurity &n_imp)
{
	m_imp = n_imp;	
	return (true);
}

double CSingleImpurityProvider::getValue(CxPosition &pos) const
{
	
	return (m_imp.getPotential(pos));
	
}

