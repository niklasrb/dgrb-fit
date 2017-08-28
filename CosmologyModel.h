#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "Constants.h"
#include "TROOT.h"
#include "TF1.h"
#include "InterpolationWrapper.h"

#include <cmath>
#include <cmath>
#include <memory>

class CosmologyModel
{
public:
	double O_m = 0; 	// matter density
	double O_dm = 0; 			// dark matter density
	double CriticalDensity = 0; // critical density
	double d_H = 0; 		// Hubble distance		
	double O_r = 0.; 			// radiation density
	double O_l = 0;			// dark energy density
	double O_k = 0;				// curvature

	
	virtual double LuminosityDistance(const double z) = 0;
	virtual double ComovingVolumeElement(const double z) = 0;
	virtual double ComovingDistance(const double z) = 0;
	virtual double HubbleRate(const double z) = 0;
};

class LambdaCDM : public CosmologyModel
{
	
private:
	std::shared_ptr<gsl1DInterpolationWrapper> ComovingDistanceSpline;
	
	const double G = 6.67e-20; // Gravtity constant in km^3 / kg*s^2
	
public:
	
	LambdaCDM(Bounds zBounds, unsigned int zGridLen)
	{
		O_m = 0.315;  			// matter density
		CriticalDensity = 3.*H_0*H_0/(8.*M_PI*G)  *(1.551812e-11) ;	// in  M_solar /  MPc^3
															// converted kg / (Mpc^2  km) to solar mass / (Mpc^3)
		d_H = c_0/H_0; 		// Hubble distance	in Mpc
		O_r = 0.; 			// radiation density
		O_l = 0.73;			// dark energy density
		O_k = 0;			// curvature
		O_dm = 0.27;
		
		std::vector<double> zGrid; zGrid.resize(zGridLen);
		std::vector<double> CD; CD.resize(zGrid.size());
		std::shared_ptr<TF1> EInverse = std::make_shared<TF1>("1/E(z) Integrand", 
											[this](double* args, double* params)  // args[0]: z
											{ return 1./sqrt( O_r*powf(1+args[0], 4) + O_m*powf(1+args[0],3) + O_k*powf(1+args[0], 2) + O_l); },
											0, 10, 0);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			zGrid.at(i) = exp( log(zBounds.first) + i*(log(zBounds.second) - log(zBounds.first))/(zGridLen-1.));
			CD.at(i) = d_H*EInverse->Integral(0, zGrid.at(i), 1e-4);
		}
		ComovingDistanceSpline = std::make_shared<gsl1DInterpolationWrapper>(zGrid.data(), zGrid.size(), CD.data(), gsl_interp_linear, 0);
		//ComovingDistanceSpline->print();
	}
	
	double ComovingVolumeElement(const double z) override
	{
		return c_0 * pow(LuminosityDistance(z),2)  /(HubbleRate(z) *pow(1.+z, 2)) ;
	}
	
	double HubbleRate(const double z)
	{
		return H_0;  // check!
	}
	
	double ComovingDistance(const double z) override
	{
		return ComovingDistanceSpline->Eval(z);
	}
	
	double LuminosityDistance(const double z) override
	{
		return ComovingDistance(z)*(1+z); 	// Only implemented for curvature =0  atm!
	}
};

#endif
