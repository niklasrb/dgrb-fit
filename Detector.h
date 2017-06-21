#ifndef DETECTOR_H
#define DETECTOR_H

#include "Constants.h"

class Detector
{
public:
	virtual double S_t_1GeV() = 0;
	virtual double DetectionEfficiency(const double S) = 0;
};

class FermiLAT : public Detector
{
	double DetectionEfficiency(const double S) override
	{
		if(S >= 5e-8_photonspercm2s1)
			return 1;
		return (S / 5e-8_photonspercm2s1)*(S / 5e-8_photonspercm2s1);
	}
	
	double S_t_1GeV()
	{
		return 5e-8_photonspercm2s1;
	}
};

#endif
