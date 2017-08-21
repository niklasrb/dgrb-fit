#ifndef DETECTOR_H
#define DETECTOR_H

#include <cmath>
#include <vector>


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
	FermiLAT() : Detector(4.5e-10) 
	{
		S_t_1 = {1e-10, 2e-10, 3e-10, 4e-10, 5e-10, 6e-10, 7e-10, 8e-10, 9e-10, 10e-10};
	}
	
	double DetectionEfficiency(const double S) override
	{
		if(S >= S_t_1GeV)
		return 1.;
	return pow(S / S_t_1GeV, 2.);
	}
};




#endif
