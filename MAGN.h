#ifndef MAGN_H
#define MAGN_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>

class MAGN : public AstrophysicalSource
{
protected:
	const double E_0 = 1._GeV; 
public:
	MAGN(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : AstrophysicalSource(_CM, tau, std::string("MAGN"))
	{
		zBounds.first = 0; zBounds.second = 6;
		GammaBounds.first = 2.37 - 2*0.32; GammaBounds.second = 2.37 + 2*0.32;
		LumBounds.first = 1e40_ergpers; LumBounds.second = 1e52_ergpers;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return pow(1+z, 2-Gamma);
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{	
		return pow(E/E_0, -Gamma);
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1./sqrt(2.*M_PI*0.32*0.32) * exp(-pow(Gamma - 2.37, 2.)/(2.*0.32*0.32));
	}			
	
	double LuminosityFunction(const double L, const double z) override
	{
		/*const double k = 3.05;*/
		const double L_151 =  pow(151./1400.,0.2) * pow(10.,-22.2/0.77) * pow(0.2*L/1._ergpers,(1./0.77));
		std::cout << "MAGN: L_151 = " << L_151 << std::endl;
		return /*k*/eta(z)*radioLuminosityFunction( L_151, z)/(0.77*1.008);
	}
	
protected:
	double eta(const double z)
	{
		return pow(c_0, 3)*pow(z*(2.+z), 2.)/(4.*pow(50.*(1. + z), 3))/CM->ComovingVolumeElement(z);
	}
	double radioLuminosityFunction(const double L, const double z)	// in erg/s
	{/*
		const double rho_lo = pow(10,-7.523);
		const double L_ls = pow(10,26.48);		// in erg/s
		const double alpha_l = 0.586;
		const double z_lo = 0.71;
		const double k_1 = 3.48;
		const double rho_ho = pow(10,-6.757);
		const double alpha_h = 2.42;
		const double L_hs = pow(10,27.39);		// in erg/s
		const double z_ho = 2.03;
		const double z_h1 = 0.568;
		const double z_h2 = 0.956;
		
		const double rho_l = rho_lo*pow((Luminosity/L_ls),-alpha_l)*exp(-Luminosity/L_ls)*pow(1.+ std::min(z, z_lo),k_1);
		if(Luminosity <= 1e12) return rho_l;
		const double rho_h = rho_ho*pow((Luminosity/L_hs),-alpha_h)*exp(-L_hs/Luminosity)*exp(-0.5*pow((z-z_ho)/(z >= z_ho ? z_h2 : z_h1), 2.));
		return rho_l + rho_h;*/
		double rho_lo = pow(10,-7.523);
		double L_ls = pow(10,26.48);
		double alpha_l = 0.586;
		double z_lo = 0.71;
		double k_1 = 3.48;
		double rho_ho = pow(10,-6.757);
		double alpha_h = 2.42;
		double L_hs = pow(10,27.39);
		double z_ho = 2.03;
		double z_h1 = 0.568;
		double z_h2 = 0.956;
		//definition of the functional dependence
		double rho_l;
		if(z<z_lo){
		rho_l = rho_lo*powf((L/L_ls),-alpha_l)*exp(-L/L_ls)*powf((1+z),k_1);
		}
		else if(z>=z_lo){
		rho_l = rho_lo*powf((L/L_ls),-alpha_l)*exp(-L/L_ls)*powf((1+z_lo),k_1);
		}
		double fh;
		if(z<z_ho){
		fh = exp(-0.5*powf(((z-z_ho)/z_h1),2));
		}
		else if(z>=z_ho){
		fh = exp(-0.5*powf(((z-z_ho)/z_h2),2));
		}
		double rho_h = rho_ho*powf((L/L_hs),-alpha_h)*exp(-L_hs/L)*fh;
		double retval = rho_l;
		if(L>1e+12){retval+=rho_h;}
		return retval;
	}
};



#endif
