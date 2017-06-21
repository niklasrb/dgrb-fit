#ifndef BLLAC_H
#define BLLAC_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>

class BLLAC : public AstrophysicalSource
{
protected:
	double E_0 = 1._GeV; // <- redefine
	double E_cut;
	double A;
	double gamma_1;
	double L_s;
	double gamma_2;
	double z_c_s;
	double alpha;
	double p_1;
	double p_2;
	
	BLLAC(std::shared_ptr<CosmologyModel> _CM, std::string name) : AstrophysicalSource(_CM, name)
	{
		m_zBounds.first = 1; m_zBounds.second = 6;
		m_GammaBounds.first = 2.1 - 2*0.26;  m_GammaBounds.second = 2.1 + 2*0.26; // check
	}
	
public:
	BLLAC(std::shared_ptr<CosmologyModel> _CM) : BLLAC(_CM, std::string("BLLAC"))
	{
		E_cut = 1; // find out
		A = 3.39e+4;
		gamma_1 = 0.27;
		L_s = 0.28e-3;
		gamma_2 = 1.86;
		z_c_s = 1.36;
		alpha = 4.53;
		p_1 = 2.24;
		p_2 = -7.37;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return powf(1 + z, 2-Gamma) * exp(-E*z/E_cut);
	}
	
	double EnergySpectrum(const double E, const double z, const double Gamma) override
	{
		return powf(E/E_0, -Gamma) * exp(-E/E_cut);
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrtf(2.*M_PI*0.26*0.26) * exp(-powf((Gamma - 2.1),2)/(2.*0.26*0.26));
	}			
	
	double LuminosityFunction(const double Luminosity, const double z, const double Gamma) override
	{
		double z_c = z_c_s*powf(Luminosity,alpha);
		double evo = 1/(powf((1+z)/(1+z_c),p_1)+powf((1+z)/(1+z_c),p_2));
		double phi_bare = A/log(10)/Luminosity/(powf(Luminosity/L_s,gamma_1)+powf(Luminosity/L_s,gamma_2));
		return phi_bare*evo;
	}
};

class LISP : public BLLAC
{
	LISP(std::shared_ptr<CosmologyModel> _CM) : BLLAC(_CM, std::string("LISP"))
	{
		E_cut = 37._GeV;
		A = 4.37e+4;
		gamma_1 = 1.19;
		L_s = 0.0308;
		gamma_2 = 0.67;
		z_c_s = 1.66;
		alpha = 0.36;
		p_1 = 4.4;
		p_2 = -2.9;
		m_GammaBounds.first = 2.08 - 2*0.15;  m_GammaBounds.second = 2.08 + 2*0.15;
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrtf(2.*M_PI*0.15*0.15) * exp(-powf((Gamma - 2.08),2)/(2.*0.15*0.15));
	}
};
	
	
	
class HSP : public BLLAC
{
	HSP(std::shared_ptr<CosmologyModel> _CM) : BLLAC(_CM, std::string("HSP"))
	{
		E_cut = 910._GeV;
		A = 98e+4;
		gamma_1 = 2.88;
		L_s = 3.15e-3;
		gamma_2 = 0.52;
		z_c_s = 4.1;
		alpha = 0.25;
		p_1 = -1.64;
		p_2 = 4.8;
		m_GammaBounds.first = 1.86 - 2*0.16;  m_GammaBounds.second = 1.86 + 2*0.16;
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrtf(2.*M_PI*0.16*0.16) * exp(-powf((Gamma - 1.86),2)/(2.*0.16*0.16));
	}		
};
#endif

