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



class AstrophysicalSource : public DGRBSource
{
	friend class Benchmark;
protected:
	Bounds zBounds;					// The redshift bounds in which the description is valid
	Bounds GammaBounds;				// The Gamma distribution in which the description is most relevant
	
	//std::function<double(const double S, const double Gamma, const double* other)> dNoverdS;	// Here an interpolation or otherwise calculated dN/dS can be stored at runtime
	
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
	
	//std::vector<std::pair<Bounds, std::shared_ptr<gsl1DInterpolationWrapper> > > APS;		// Contains pair of Energy bin and C_p(S_t_1) Spline
	std::shared_ptr<AstrophysicalSourceAPS<std::shared_ptr<gsl1DInterpolationWrapper> > > APS;
	
	void printResults(double S_t_1);
	
	TGraph MakeGraph()
	{
		std::vector<double> x, y;
		for(unsigned int i = 0; i < Intensity.size(); i++) 
		{
			x.push_back((Intensity[i].first.first + Intensity[i].first.second)/2.);
			y.push_back(Intensity[i].second);
			//std::cout << x[i] << '\t' << y[i] << std::endl;
		}
		return TGraph(x.size(), x.data(), y.data());
	}
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

void AstrophysicalSource::printResults(double S_t_1)
{
	std::cout << Name <<  "  - Bin: Intensity" << std::endl;
	for(unsigned int i = 0; i < Intensity.size(); i++) 
		std::cout << "[" << Intensity.at(i).first.first << "," << Intensity.at(i).first.second <<"]:" << Intensity.at(i).second << '\t';
	std::cout << std::endl;// << "Bin: C_p  ";
	//for(unsigned int i = 0; i < APS.size(); i++) 
	//	std::cout << "[" << APS.at(i).first.first << "," << APS.at(i).first.second <<"]:" << APS.at(i).second->Eval(S_t_1) << '\t';
	//std::cout << std::endl;
}

#include "MAGN.h"
#include "FSRQ.h"
#include "BLLAC.h"
#include "SFG.h"

#endif
