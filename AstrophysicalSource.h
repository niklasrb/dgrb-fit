#ifndef ASTROPHYSICALSOURCE_H
#define ASTROPHYSICALSOURCE_H

/// This header contains a model class for astrophysical sources
/// Basically sources where luminosity plays the most important role
/// Most of the functions are meant to be overriden for specific sources


#include <functional>
#include <cmath>
#include <memory>

#include "Constants.h"
#include "Source.h"




class AstrophysicalSource : public DGRBSource
{
	friend class Benchmark;
protected:
	Bounds zBounds;					// The redshift bounds in which the description is valid
	Bounds GammaBounds;				// The Gamma distribution in which the description is most relevant
	
	std::function<double(const double S, const double Gamma, const double* other)> dNoverdS;	// Here an interpolation or otherwise calculated dN/dS can be stored at runtime
	
	AstrophysicalSource(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : DGRBSource(CM, tau) {}
	AstrophysicalSource(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string name) : DGRBSource(CM, tau, name) {}
	~AstrophysicalSource() {}

public:
	// Functions that have to be modelled for the different source classes
	virtual double kCorrection(const double E, const double z, const double Gamma) = 0;
	virtual double literal_F(const double E, const double z, const double Gamma) = 0;
	virtual double GammaDistribution(const double Gamma) = 0;											// dN/dGamma
	virtual double LuminosityFunction(const double Luminosity, const double z, const double Gamma)= 0;
	
	// These functions are already implemented for convenience, but can be overriden
	virtual double EnergySpectrum(const double E, const double z, const double Gamma);					// dN/dE
	virtual double RescaledLuminosityFunction(const double Luminosity, const double z, const double Gamma);
	
	
};

// Preimplemented functions that use the other implementations
double AstrophysicalSource::RescaledLuminosityFunction(const double Luminosity, const double z, const double Gamma)
{
	return LuminosityFunction(Luminosity, z, Gamma)*GammaDistribution( Gamma);
}

double AstrophysicalSource::EnergySpectrum(const double E, const double z, const double Gamma)
{
	return literal_F(E, z, Gamma)*exp(-(*tau)(E, z));
}

#include "MAGN.h"
#include "FSRQ.h"
#include "BLLAC.h"
#include "SFG.h"

#endif
