#ifndef CSINGLEIMPURITYPROVIDER_H_
#define CSINGLEIMPURITYPROVIDER_H_

#include <utils/IValueProvider.h>
#include <geometry/CxImpurity.hpp>

class CSingleImpurityProvider : public IValueProvider
{
public:
	virtual ~CSingleImpurityProvider();
	static CSingleImpurityProvider * getInstance(){return (m_instance);};
	virtual double getValue(CxPosition &) const;
	bool setValue(CxImpurity &n_imp);
protected:
	CSingleImpurityProvider();
	static CSingleImpurityProvider * m_instance;
	CxImpurity m_imp;
};

#endif /*CSINGLEIMPURITYPROVIDER_H_*/
