#ifndef SOURCE_H
#define SOURCE_H

/// This header contains the model class for all possible DGRB contributions
/// The functions Intensity and Autocorrelation are meant to be calculated during runtime
/// This can be done by different algorithms, this is just the class where the result is supposed to be saved

#include <iostream>
#include <functional>
#include <utility>
#include <memory>

#include "Constants.h"
#include "CosmologyModel.h"
#include "EBLAbsorbtionCoefficient.h"
#include "AngularPowerSpectrum.h"


/* This class is the base class for all source populations
 * Containing a window function and a source density FT
 * the EBL absorption coefficient and Cosmology Model
 * Also holds the intensity values
 */
class DGRBSource
{
public:
	std::string Name;
	
	DGRBSource(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : DGRBSource(_CM, tau, std::string("")) {} 				// simple constructors
	DGRBSource(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string _name) : Name(_name), CM(_CM), tau(tau) {}	// standard copy constructors should work
	
	std::shared_ptr<CosmologyModel> CM;							// A pointer to a class containing the parameters of a cosmology model
	std::shared_ptr<EBLAbsorbtionCoefficient> tau;				// A pointer to a class modelling the absorbtion due to Extragalactic Background Light, commonly tau
	
	
	// These functions have to be set at runtime
	std::function<double(const double E, const double z)> WindowFunction = 0;					// All sources have window functions
	std::function<double(const double k, const double M, const double z)> SourceDensityFT = 0;	// The fourier transform of the Source Density, important for calculating cross sections with Galaxy Catalogues
	
	// The calculated results will be saved here
	std::vector<double > Intensity;		// A list of  Intensities, depending on EnergyBin  (maybe encapsulate this better)

};



#endif
