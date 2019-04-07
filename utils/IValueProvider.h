#ifndef IVALUEPROVIDER_H_
#define IVALUEPROVIDER_H_
#include <geometry/CxPosition.hpp>

class IValueProvider
{
public:
	virtual ~IValueProvider();
	virtual double getValue(CxPosition &) const = 0;
};

#endif /*IVALUEPROVIDER_H_*/
