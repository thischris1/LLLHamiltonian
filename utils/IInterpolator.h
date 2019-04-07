#ifndef IINTERPOLATOR_H_
#define IINTERPOLATOR_H_
#include <geometry/CxPosition.hpp>
class IInterpolator
{
public:
	virtual double getValue(const CxPosition &zielPunkt, std::vector<CxPosition> Punkte, std::vector<double> values) = 0;
	virtual ~IInterpolator();
};

#endif /*IINTERPOLATOR_H_*/
