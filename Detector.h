#ifndef DETECTOR_H
#define DETECTOR_H

#include <cmath>
#include <vector>
#include "Constants.h"

/* This is a small class that contains the detection efficiencies
 * and the flux threshold > 1GeV
 * Since the flux threshold depends more on the galaxy catalog
 * this class should be abolished in favor of the galaxy catalog class
 */
class Detector
{
public:
	const double S_t_1GeV;
	virtual double DetectionEfficiency(const double S) {return 0;}
	Detector(double S_t_1GeV) : S_t_1GeV(S_t_1GeV) {}
	~Detector() {}
	std::vector<double> S_t_1;
};

class FermiLAT : public Detector
{
public:
	FermiLAT() : Detector(5e-10_photonspercm2s) 
	{
		S_t_1 = {1e-10_photonspercm2s, 2e-10_photonspercm2s, 3e-10_photonspercm2s, 4e-10_photonspercm2s, 
					5e-10_photonspercm2s, 6e-10_photonspercm2s, 7e-10_photonspercm2s, 8e-10_photonspercm2s, 9e-10_photonspercm2s, 10e-10_photonspercm2s};
	}
	
	double DetectionEfficiency(const double S) override
	{
		if(S >= 5e-8_photonspercm2s)
			return 1.;
		return pow(S / 5e-8_photonspercm2s, 2.);
	}
};




#endif
