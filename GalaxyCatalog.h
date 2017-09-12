#ifndef GALAXYCATALOG_H
#define GALAXYCATALOG

#include <cmath>
#include <functional>
#include <string>

/* These classes model Galaxy catalogues
 * They are not used at the moment, but were intended for use in cross correlation computations
 * Maybe make them take on detection efficiency
 * But those are used differently, so I'm not sure
 */


class GalaxyCatalog
{
protected:
	double M_th = 1;
	double sigma_logM = 1;
	double M_cut = 1;
	double M_1 = 1;
	double alpha = 0;
	double beta = 1;
	double m = 0;
	double z_0 = 1;
	//double S_t_1 = 0;
	
public:
	std::string name = "";
	
	double N_central(const double M);
	double N_satellite(const double M);
	
	double WindowFunction(const double z);
	std::function<double(const double r, const double M, const double z)> SourceDensity = 0;
	std::function<double(const double k, const double M, const double z)> SourceDensityFT = 0;
	std::function<double(const double z)> EffectiveGalaxyDensity = 0;
	
	/*double S_t_1GeV() { return S_t_1; }
	
	double DetectionEfficiency(const double S)
	{
		return DetectionEfficiency(S, S_t_1);
	}
	
	static double DetectionEfficiency(const double S, const double S_t_1Gev);*/
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

/*double GalaxyCatalog::DetectionEfficiency(const double S, const double S_t_1Gev)
{
	if(S >= S_t_1Gev)
		return 1;
	return pow(S / S_t_1Gev, 2.);
}*/



#endif
