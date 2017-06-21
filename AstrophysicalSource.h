#ifndef ASTROPHYSICALSOURCE_H
#define ASTROPHYSICALSOURCE_H


#include <functional>
#include <cmath>

#include "Constants.h"
#include "Source.h"




class AstrophysicalSource : public DGRBSource
{
	friend class Benchmark;
protected:
	Bounds m_zBounds;
	Bounds m_GammaBounds;
	std::function<double(const double *)> dNoverdS;
	
public:
	

	AstrophysicalSource(std::shared_ptr<CosmologyModel> _CM) : DGRBSource(_CM) {}
	AstrophysicalSource(std::shared_ptr<CosmologyModel> _CM, std::string _name) : DGRBSource(_CM, _name) {}
	~AstrophysicalSource() {}
	
	virtual double kCorrection(const double E, const double z, const double Gamma) = 0;
	virtual double EnergySpectrum(const double E, const double z, const double Gamma)= 0;					// dN/dE
	virtual double GammaDistribution(const double Gamma) = 0;											// dN/dGamma
	virtual double LuminosityFunction(const double Luminosity, const double z, const double Gamma)= 0;
	virtual double RescaledLuminosityFunction(const double Luminosity, const double z, const double Gamma);
	
	
};

double AstrophysicalSource::RescaledLuminosityFunction(const double Luminosity, const double z, const double Gamma)
{
	return LuminosityFunction(Luminosity, z, Gamma)*GammaDistribution( Gamma);
}


#include "MAGN.h"
#include "FSRQ.h"
#include "BLLAC.h"

#endif
