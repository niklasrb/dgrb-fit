#ifndef MAGN_H
#define MAGN_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>

class MAGN : public AstrophysicalSource
{
protected:
	double m_E0 = 1* GeV; // <- redefine
public:
	MAGN(std::shared_ptr<CosmologyModel> _CM) : AstrophysicalSource(_CM, std::string("MAGN"))
	{
		m_zBounds.first = 1; m_zBounds.second = 6;
		m_GammaBounds.first = 2.37 - 2*0.32; m_GammaBounds.second = 2.37 + 2*0.32;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return powf(1+z, 2-Gamma);
	}
	
	double EnergySpectrum(const double E, const double z, const double Gamma) override
	{
		double tau = 0; 		// Implement
		return powf(E/m_E0, -Gamma)*exp(-tau);
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrtf(2.*M_PI*0.32*0.32) * exp(-(Gamma - 2.37)*(Gamma - 2.37)/(2.*0.32*0.32));
	}			
	
	double LuminosityFunction(const double Luminosity, const double z, const double Gamma) override
	{
		return eta(z)*radioLuminosityFunction( powf(151./1200., 0.2)* powf(10., -22.2/0.77) * powf(0.2*Luminosity, 1./0.77), z)/(0.77*1.008);
	}
	
protected:
	double eta(const double z)
	{
		return powf(c_0, 3)*powf(z*(2+z), 2)/(4.*powf(50.*(1. + z), 3))/CM->ComovingVolume(z);
	}
	double radioLuminosityFunction(const double Luminosity, const double z)
	{
		const double rho_lo = powf(10,-7.523);
		const double L_ls = powf(10,26.48);
		const double alpha_l = 0.586;
		const double z_lo = 0.71;
		const double k_1 = 3.48;
		const double rho_ho = powf(10,-6.757);
		const double alpha_h = 2.42;
		const double L_hs = powf(10,27.39);
		const double z_ho = 2.03;
		const double z_h1 = 0.568;
		const double z_h2 = 0.956;
		
		//std::cout << powf((Luminosity/L_ls),-alpha_l) << '\t' << exp(-Luminosity/L_ls) << '\t' << powf(1+ std::min(z, z_lo),k_1) << std::endl;
		double rho_l = rho_lo*powf((Luminosity/L_ls),-alpha_l)*exp(-Luminosity/L_ls)*powf(1+ std::min(z, z_lo),k_1);
		if(Luminosity <= 1e12) return rho_l;
		double rho_h = rho_ho*powf((Luminosity/L_hs),-alpha_h)*exp(-L_hs/Luminosity)*exp(-0.5*powf((z-z_ho)/(z >= z_ho ? z_h2 : z_h1), 2));
		return rho_l + rho_h;
	}
};



#endif
