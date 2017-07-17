#ifndef GALAXYCATALOG_H
#define GALAXYCATALOG

#include <cmath>
#include <functional>
#include <string>

class GalaxyCatalog
{
protected:
	const double M_th = 1;
	const double sigma_logM = 1;
	const double M_cut = 1;
	const double M_1 = 1;
	const double alpha = 0;
	const double beta = 1;
	const double m = 0;
	const double z_0 = 1;
	
public:
	std::string name = "";
	
	double N_central(const double M);
	double N_satellite(const double M);
	
	double WindowFunction(const double z);
	std::function<double(const double r, const double M, const double z)> SourceDensity = 0;
	std::function<double(const double k, const double M, const double z)> SourceDensityFT = 0;
	std::function<double(const double z)> EffectiveGalaxyDensity = 0;


};

double GalaxyCatalog::WindowFunction(const double z)
{
	return beta / tgamma((m+1.)/beta) * pow(z, m) / pow(z_0, m+1.) * exp(-pow(z/z_0, beta));
}

double GalaxyCatalog::N_central(const double M)
{
	return (1. + erf( (log(M) - log(M_th))/sigma_logM))/2.;
}

double GalaxyCatalog::N_satellite(const double M)
{
	return powf(M/M_1, 2) * exp(-M/M_cut);
}



#endif
