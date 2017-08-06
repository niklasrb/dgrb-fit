#ifndef EBLABSORBTIONCOEFFICIENT_H
#define EBLABSORBTIONCOEFFICIENT_H

#include <cmath>
#include "InterpolationWrapper.h"


class EBLAbsorbtionCoefficient
{
protected:
	gsl2DInterpolationWrapper interp;

public:

	EBLAbsorbtionCoefficient(gsl2DInterpolationWrapper interp) : interp(interp) {}

	double operator ()(const double E, const double z)
	{
		return interp.Eval(E, z);
	}
};



#endif
