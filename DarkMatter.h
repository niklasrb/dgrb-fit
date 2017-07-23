#ifndef DARKMATTER_H
#define DARKMATTER_H

/// This header contains the classes for different Dark Matter models
/// They also inherit from 

#include <string>
#include <memory>
#include <cmath>

#include "Constants.h"
#include "Source.h"
#include "HaloModel.h"


class DarkMatter : public DGRBSource
{
protected:
	std::shared_ptr<gsl1DInterpolationWrapper> dNdLogx;

public:
	std::shared_ptr<HaloModel> HM;
	const double Mass;
	std::function<double(const double E, const double z)> EnergySpectrum = 0;  // dN/dE
	
	DarkMatter(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::shared_ptr<gsl1DInterpolationWrapper> dNdLogx, const double m, std::string name) : DGRBSource(_CM, tau, name), dNdLogx(dNdLogx), HM(HM), Mass(m)  {}

};

class AnnihilatingDM : public DarkMatter
{
public:
	const double ThermalizedAnnihilationCrossSection;
	
	AnnihilatingDM(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::shared_ptr<gsl1DInterpolationWrapper> dNdLogx, const double m, const double ThermalizedAnnihilationCrossSection) : DarkMatter(_CM, HM, tau, dNdLogx, m, "Annihilating Dark Matter"), ThermalizedAnnihilationCrossSection(ThermalizedAnnihilationCrossSection)
	{
		// Set EnergySpectrum
		EnergySpectrum = [this] ( const double E, const double z)
						{ 	//std::cout << log(2.*E/(e*c_0*c_0*Mass)) << '\t' << this->dNdLogx->Eval(log(2.*E/(e*c_0*c_0*Mass)))/E << std::endl;
							return this->dNdLogx->Eval(log(2.*E/(e*c_0*c_0*Mass))) /E; };
		
		WindowFunction = [this] ( const double E, const double z) 
						{ std::cout << this->ThermalizedAnnihilationCrossSection / (8.*M_PI)  << '\t' <<  pow(CM->O_dm * CM->CriticalDensity / Mass, 2.) * pow(1.+z, 3.) << '\t' <<  this->EnergySpectrum(E*(1.+z), z) << '\t' << exp(-(*(this->tau))(E,z)) << std::endl;
							return this->ThermalizedAnnihilationCrossSection / (8.*M_PI)  * pow(CM->O_dm * CM->CriticalDensity / Mass, 2.) * pow(1.+z, 3.) * this->EnergySpectrum(E*(1.+z), z) * exp(-(*(this->tau))(E,z)); };
	
		SourceDensityFT = [this] (const double k, const double M, const double z)
						{ return this->HM->SourceDensitySubhaloBoostFT(k, M, z) / pow(this->HM->ClumpingFactor(z), 2.); };
	}
	
	
};

class DecayingDM : public DarkMatter
{
public:
	const double HalfLife;
	
	DecayingDM(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::shared_ptr<gsl1DInterpolationWrapper> dNdLogx, const double m, const double HalfLife) : DarkMatter(_CM, HM, tau, dNdLogx, m, "Decaying Dark Matter"), HalfLife(HalfLife)
	{
		EnergySpectrum = [this] ( const double E, const double z)
						{ return this->dNdLogx->Eval(log(1.*E/(e*c_0*c_0*Mass))) /E; };
		
		WindowFunction = [this] ( const double E, const double z) 
						{  return 1. / (4.*M_PI)  * CM->O_dm * CM->CriticalDensity / (Mass*this->HalfLife)* this->EnergySpectrum(E*(1.+z), z) * exp(-(*(this->tau))(E,z)); };
	
		SourceDensityFT = HM->NFWHaloDensityProfileFT;
	}
	
};


#endif
