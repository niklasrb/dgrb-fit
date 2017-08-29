#ifndef FSRQ_H
#define FSRQ_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>

class FSRQ : public AstrophysicalSource
{
protected:
	const double E_cut = 6._GeV;
	const double E_0 = 1._GeV;
	
public:
	FSRQ(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : AstrophysicalSource(_CM, tau, std::string("FSRQ"))
	{
		zBounds.first = 0; zBounds.second = 6;
		GammaBounds.first = 2.44 - 2*0.18; GammaBounds.second = 2.44 + 2*0.18;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return pow(1 + z, 2 - Gamma) * exp(-(sqrt(E*(1+z)) - E) / sqrt(E_cut));
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{
		return pow(E/E_0, -Gamma) * exp(- sqrt(E/ E_cut));
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1./sqrt(2.*M_PI*0.18*0.18) * exp(-pow((Gamma - 2.44), 2)/(2.*0.18*0.18));
	}			
	
	double LuminosityFunction(const double L, const double z, const double Gamma) override
	{
		//defining parameters of the best fit model
		const double A = 3.06e-9/pow(1._Mpc,3)*1e48_ergpers; 
		const double gamma_1 = 0.21;
		const double L_s = 0.84e48_ergpers;
		const double gamma_2 = 1.58;
		const double z_c_s = 1.47;
		const double alpha = 0.21;
		const double p_1 = 7.35;
		const double p_2 = -6.51;
		
		double z_c = z_c_s*pow(L/1e48_ergpers, alpha);
		double evo = 1./(pow((1.+z)/(1+z_c),p_1)+pow((1.+z)/(1.+z_c),p_2));
		double phi_bare = A/(log(10.)*L*(pow(L/L_s,gamma_1)+pow(L/L_s,gamma_2)));
		return phi_bare*evo;
	}
	

};



#endif
