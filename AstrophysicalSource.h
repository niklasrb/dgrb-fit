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
#include "TGraph.h"


/* Describes the astrophysical source populations
 * With all their characteristic functions
 * The APS is saved as 1D interpolation objects
 */

class AstrophysicalSource : public DGRBSource
{
	friend class Benchmark;
protected:
	Bounds zBounds;					// The redshift bounds in which the description is valid
	Bounds GammaBounds;				// The Gamma distribution in which the description is most relevant
	Bounds LumBounds;				// luminosity bounds
	
	AstrophysicalSource(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : DGRBSource(CM, tau) {}
	AstrophysicalSource(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string name); 
	~AstrophysicalSource() {}

public:
	// Functions that have to be modelled for the different source classes
	virtual double kCorrection(const double E, const double z, const double Gamma) = 0;					// Energy 
	virtual double literal_F(const double E, const double z, const double Gamma) = 0;
	virtual double GammaDistribution(const double Gamma) = 0;											// dN/dGamma
	virtual double LuminosityFunction(const double Luminosity, const double z)= 0;						
	
	// These functions are already implemented for convenience, but can be overriden
	virtual double EnergySpectrum(const double E, const double z, const double Gamma);					// dN/dE
	virtual double RescaledLuminosityFunction(const double Luminosity, const double z, const double Gamma);
	virtual double EnergySpectrumOverK(const double E, const double z, const double Gamma);				// dN/dE / k correction
	

	std::shared_ptr<AstrophysicalSourceAPS<std::shared_ptr<gsl1DInterpolationWrapper> > > APS;			// 2dim matrix that will hold the APS 1D interpolation objects
	
	
};

AstrophysicalSource::AstrophysicalSource(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string name)
:	DGRBSource(CM, tau, name)
{
	zBounds = Bounds(0, 6);
	GammaBounds = Bounds(0., 0.);
	LumBounds = Bounds(1e40_ergpers, 1e52_ergpers);
}

// Preimplemented functions that use the other implementations
double AstrophysicalSource::RescaledLuminosityFunction(const double Luminosity, const double z, const double Gamma)
{
	return LuminosityFunction(Luminosity, z)*GammaDistribution( Gamma);
}

double AstrophysicalSource::EnergySpectrum(const double E, const double z, const double Gamma)
{
	return literal_F(E, z, Gamma)*exp(-(*tau)(E, z));
}

double AstrophysicalSource::EnergySpectrumOverK(const double E, const double z, const double Gamma)
{
	return EnergySpectrum(E, z, Gamma)/kCorrection(E, z, Gamma);
}



#include "MAGN.h"
#include "FSRQ.h"
#include "BLLAC.h"
#include "SFG.h"

#endif
