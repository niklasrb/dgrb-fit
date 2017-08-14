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
	MAGN(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : AstrophysicalSource(_CM, tau, std::string("MAGN"))
	{
		zBounds.first = 1; zBounds.second = 6;
		GammaBounds.first = 2.37 - 2*0.32; GammaBounds.second = 2.37 + 2*0.32;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return pow(1+z, 2-Gamma);
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{
		double tau = 0; 		// Implement
		return pow(E/m_E0, -Gamma)*exp(-tau);
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrt(2.*M_PI*0.32*0.32) * exp(-(Gamma - 2.37)*(Gamma - 2.37)/(2.*0.32*0.32));
	}			
	
	double LuminosityFunction(const double Luminosity, const double z, const double Gamma) override
	{
		return eta(z)*radioLuminosityFunction( pow(151./1200., 0.2)* pow(10., -22.2/0.77) * pow(0.2*Luminosity, 1./0.77), z)/(0.77*1.008);
	}
	
protected:
	double eta(const double z)
	{
		return pow(c_0, 3)*pow(z*(2+z), 2)/(4.*pow(50.*(1. + z), 3))/CM->ComovingVolume(z);
	}
	double radioLuminosityFunction(const double Luminosity, const double z)
	{
		const double rho_lo = pow(10,-7.523);
		const double L_ls = pow(10,26.48);
		const double alpha_l = 0.586;
		const double z_lo = 0.71;
		const double k_1 = 3.48;
		const double rho_ho = pow(10,-6.757);
		const double alpha_h = 2.42;
		const double L_hs = pow(10,27.39);
		const double z_ho = 2.03;
		const double z_h1 = 0.568;
		const double z_h2 = 0.956;
		
		//std::cout << powf((Luminosity/L_ls),-alpha_l) << '\t' << exp(-Luminosity/L_ls) << '\t' << powf(1+ std::min(z, z_lo),k_1) << std::endl;
		double rho_l = rho_lo*pow((Luminosity/L_ls),-alpha_l)*exp(-Luminosity/L_ls)*pow(1+ std::min(z, z_lo),k_1);
		if(Luminosity <= 1e12) return rho_l;
		double rho_h = rho_ho*pow((Luminosity/L_hs),-alpha_h)*exp(-L_hs/Luminosity)*exp(-0.5*pow((z-z_ho)/(z >= z_ho ? z_h2 : z_h1), 2.));
		return rho_l + rho_h;
	}
};



#endif
