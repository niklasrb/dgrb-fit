#ifndef LINEARMATTERPOWERSPECTRUM_H
#define LINEARMATTERPOWERSPECTRUM_H


#include "gsl2DInterpolationWrapper.h"


class LinearMatterPowerSpectrum
{
protected:
	gsl2DInterpolationWrapper interp;
	
public:
	LinearMatterPowerSpectrum(gsl2DInterpolationWrapper interp) : interp(interp) {}
	
	double operator ()(const double k, const double z)
	{
		return interp.Eval(k, z);
	}
};


#endif
