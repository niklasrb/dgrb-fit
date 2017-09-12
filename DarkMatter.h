#ifndef DARKMATTER_H
#define DARKMATTER_H

/// This header contains the classes for different Dark Matter models
/// They also inherit from 

#include <string>
#include <memory>
#include <cmath>
#include <vector>

#include "Constants.h"
#include "Source.h"
#include "HaloModel.h"
#include "AngularPowerSpectrum.h"

/* Classes that hold information about Dark Matter objects
 * Energy spectrum is taken from simulations
 * Parameter here is Mass
 * APS is just a 3D matrix pointing to double values
 */
class DarkMatter : public DGRBSource
{
protected:
	std::shared_ptr<gsl2DInterpolationWrapper> dNdLogx;

public:
	std::shared_ptr<HaloModel> HM;
	const double Mass;
	
	std::shared_ptr<AngularPowerSpectrum<double> > APS;
	
	std::function<double(const double E, const double z)> EnergySpectrum = 0;  // dN/dE
	
	DarkMatter(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::shared_ptr<gsl2DInterpolationWrapper> dNdLogx, const double m, std::string name) : DGRBSource(_CM, tau, name), dNdLogx(dNdLogx), HM(HM), Mass(m)  {}

};

/// Annihilating Dark matter is characterized by the Thermalized Annihilation Cross Section
class AnnihilatingDM : public DarkMatter
{
public:
	const double ThermalizedAnnihilationCrossSection;
	
	AnnihilatingDM(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::shared_ptr<gsl2DInterpolationWrapper> dNdLogx, const double m, const double ThermalizedAnnihilationCrossSection) : DarkMatter(_CM, HM, tau, dNdLogx, m, "Annihilating Dark Matter"), ThermalizedAnnihilationCrossSection(ThermalizedAnnihilationCrossSection)
	{
		// Set EnergySpectrum
		EnergySpectrum = [dNdLogx, m] ( const double E, const double z)
						{ 	return std::max(0., dNdLogx->Eval(m, log10(2.*E/m))) /E; };
		
		WindowFunction = [this] ( const double E, const double z) 
						{ 	return this->ThermalizedAnnihilationCrossSection / (8.*M_PI)  * pow(CM->O_dm * CM->CriticalDensity / Mass, 2.) * pow(1.+z, 3.) * this->HM->ClumpingFactor(z) * this->EnergySpectrum(E*(1.+z), z) * exp(-(*(this->tau))(E*(1.+z),z)); };
	
		SourceDensityFT = [this] (const double k, const double M, const double z)	// Load the source density from the Halo Model		// This should probably be moved out of the HaloModel class, and computed internally in the DM class
						{ std::cout << "SDSBFT: " << this->HM->SourceDensitySubhaloBoostFT(k, M, z) << " CF: " << this->HM->ClumpingFactor(z) << std::endl;
							return this->HM->SourceDensitySubhaloBoostFT(k, M, z) / this->HM->ClumpingFactor(z); };
	}						
	
	
};

/// Decaying DM is characterized by its Half Life
class DecayingDM : public DarkMatter
{
public:
	const double HalfLife;
	
	DecayingDM(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::shared_ptr<gsl2DInterpolationWrapper> dNdLogx, const double m, const double HalfLife) : DarkMatter(_CM, HM, tau, dNdLogx, m, "Decaying Dark Matter"), HalfLife(HalfLife)
	{
		EnergySpectrum = [dNdLogx, m] ( const double E, const double z)
						{ return std::max(0., dNdLogx->Eval(m, log10(E/m))) /E; };
		
		WindowFunction = [this] ( const double E, const double z) 
						{  return 1. / (4.*M_PI)  * CM->O_dm * CM->CriticalDensity / (Mass*this->HalfLife)* this->EnergySpectrum(E*(1.+z), z) * exp(-(*(this->tau))(E*(1.+z),z)); };
	
		SourceDensityFT = HM->NFWHaloDensityProfileFT;
	}
	
};


#endif
