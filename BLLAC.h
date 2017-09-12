#ifndef BLLAC_H
#define BLLAC_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>
#include <memory>

/* These classes model Blazar contributions
 * Low and Intermediate synchroton peaked LISP
 * and high synchroton peaked HSP 
 */

class BLLAC : public AstrophysicalSource
{
protected:
	double E_0 = 1._GeV;
	double E_cut;
	double A;
	double gamma_1;
	double L_s;
	double gamma_2;
	double z_c_s;
	double alpha;
	double p_1;
	double p_2;

	BLLAC(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string name, const double E_cut) : AstrophysicalSource(CM, tau, name), E_cut(E_cut)
	{
		zBounds.first = 0; zBounds.second = 6;
		GammaBounds.first = 2.1 - 2*0.26;  GammaBounds.second = 2.1 + 2*0.26; 
	}
	
public:
	BLLAC(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, const double E_cut) : BLLAC(CM, tau, std::string("BLLAC"), E_cut)
	{
		A = 3.39e-9/pow(1._Mpc,3)*1e48_ergpers;
		gamma_1 = 0.27;
		L_s = 0.28*1e45_ergpers;
		gamma_2 = 1.86;
		z_c_s = 1.36;
		alpha = 4.53;
		p_1 = 2.24;
		p_2 = -7.37;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		return pow(1 + z, 2-Gamma) * exp(-E*z/E_cut);
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{
		return pow(E/E_0, -Gamma) * exp(-E/E_cut);
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrt(2.*M_PI*0.26*0.26) * exp(-pow((Gamma - 2.1),2)/(2.*0.26*0.26));
	}			
	
	double LuminosityFunction(const double L, const double z) override
	{
		const double z_c = z_c_s*pow(L/1e48_ergpers,alpha);
		const double evo = 1./(pow((1.+z)/(1.+z_c),p_1) + pow((1+z)/(1+z_c),p_2));
		const double phi_bare = A/log(10)/L/(pow(L/L_s,gamma_1)+pow(L/L_s,gamma_2));
		return phi_bare*evo;
	}
};

class LISP : public BLLAC
{
public:
	LISP(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, const double E_cut = 37._GeV) : BLLAC(CM, tau, std::string("LISP"), E_cut)
	{
		A = 4.37e-9/pow(1._Mpc,3)*1e48_ergpers;
		gamma_1 = 1.19;
		L_s = 30.8e45_ergpers;
		gamma_2 = 0.67;
		z_c_s = 1.66;
		alpha = 0.36;
		p_1 = 4.4;
		p_2 = -2.9;
		GammaBounds.first = 2.08 - 2*0.15;  GammaBounds.second = 2.08 + 2*0.15;
		zBounds.first = 0; zBounds.second = 3.;
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrt(2.*M_PI*0.15*0.15) * exp(-pow((Gamma - 2.08),2)/(2.*0.15*0.15));
	}
};
	
	
	
class HSP : public BLLAC
{
public:
	HSP(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, const double E_cut = 910._GeV) : BLLAC(CM, tau, std::string("HSP"), E_cut)
	{
		A = 98.e-9/pow(1._Mpc,3)*1e48_ergpers;
		gamma_1 = 2.88;
		L_s = 3.15e45_ergpers;
		gamma_2 = 0.52;
		z_c_s = 4.1;
		alpha = 0.25;
		p_1 = -1.64;
		p_2 = 4.8;
		GammaBounds.first = 1.86 - 2*0.16;  GammaBounds.second = 1.86 + 2*0.16;
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1/sqrt(2.*M_PI*0.16*0.16) * exp(-pow((Gamma - 1.86),2)/(2.*0.16*0.16));
	}		
};
#endif

