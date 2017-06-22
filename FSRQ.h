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
	FSRQ(std::shared_ptr<CosmologyModel> _CM) : AstrophysicalSource(_CM, std::string("FSRQ"))
	{
		m_zBounds.first = 1; m_zBounds.second = 6;
		m_GammaBounds.first = 2.44 - 2*0.18; m_GammaBounds.second = 2.44 + 2*0.18;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return powf(1 + z, 2 - Gamma) * exp(-(sqrt(E*(1+z)) - E) / sqrt(E_cut));
	}
	
	double EnergySpectrum(const double E, const double z, const double Gamma) override
	{
		return powf(E/E_0, -Gamma) * exp(- sqrt(E/ E_cut));
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrtf(2.*M_PI*0.18*0.18) * exp(-powf((Gamma - 2.44), 2)/(2.*0.18*0.18));
	}			
	
	double LuminosityFunction(const double Luminosity, const double z, const double Gamma) override
	{
		//defining parameters of the best fit model
		const double A = 3.06e+4; 
		const double gamma_1 = 0.21;
		const double L_s = 0.84;
		const double gamma_2 = 1.58;
		const double z_c_s = 1.47;
		const double alpha = 0.21;
		const double p_1 = 7.35;
		const double p_2 = -6.51;
		
		double z_c = z_c_s*powf(Luminosity, alpha);
		double evo = 1/(powf((1+z)/(1+z_c),p_1)+powf((1+z)/(1+z_c),p_2));
		double phi_bare = A/log(10)/Luminosity/(powf(Luminosity/L_s,gamma_1)+powf(Luminosity/L_s,gamma_2));
		return phi_bare*evo;
	}
	

};



#endif
