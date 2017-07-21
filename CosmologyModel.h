#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "Constants.h"
#include <cmath>
#include "TROOT.h"
#include "TF1.h"
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
	virtual double ComovingVolume(const double z) = 0;
	virtual double ComovingDistance(const double z) = 0;
	virtual double HubbleRate(const double z) = 0;
};

class LambdaCDM : public CosmologyModel
{
	
private:
	std::shared_ptr<TF1> EInverse = std::make_shared<TF1>("1/E(z) Integrand", 
											[this](double* args, double* params)  // args[0]: z
											{ return 1./sqrt( O_r*powf(1+args[0], 4) + O_m*powf(1+args[0],3) + O_k*powf(1+args[0], 2) + O_l); },
											0, 10, 0);
	
public:
	
	
	LambdaCDM()
	{
		O_m = 0.27;  			// matter density
		CriticalDensity = 3*H_0*H_0/(8*M_PI*G);
		d_H = c_0/H_0; 		// Hubble distance		
		O_r = 0.; 			// radiation density
		O_l = 0.73;			// dark energy density
		O_k = 0;			// curvature
		O_dm = 0.23;
	}
	
	double ComovingVolume(const double z) override
	{
		return c_0 * LuminosityDistance(z)  /(H_0 *(1.+z)*(1.+z)) ;
		//double coeff = -0.0002+1.013*z+0.742*powf(z,2)-0.217*powf(z,3)+0.0343*powf(z,4)-0.00211*powf(z,5);
		//return 4*M_PI*powf(c_0/H_0,3)*powf(coeff,2)/powf((1+z),2)/sqrt(O_l+O_m*powf((1+z),3));
	}
	
	double HubbleRate(const double z)
	{
		return H_0;  // check!
	}
	
	double ComovingDistance(const double z) override
	{
		return d_H*EInverse->Integral(0, z, 1e-3);
	}
	
	double LuminosityDistance(const double z) override
	{
		return ComovingDistance(z)*(1+z); 	// Only implemented for curvature =0  atm!
	}
};

#endif
