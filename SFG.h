#ifndef SFG_H
#define SFG_H

#include "AstrophysicalSource.h"
#include "Constants.h"

#include <cmath>
#include <algorithm>

class SFG : public AstrophysicalSource
{
protected:
	double Gamma_X;
	double alpha;
	double beta;
	double sigma;
	double k_L1;
	double k_L2;
	double z_L;
	double k_Phi1;
	double k_Phi2;
	double z_Phi;
	double L_IR_0; // maybe obsolete
	double Phi_s_0;
	
	SFG(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string name) : AstrophysicalSource(_CM, tau, name)
	{
		zBounds.first = 1; zBounds.second = 6;
		GammaBounds.first = 0; GammaBounds.second = 0;
	}
	
public:
	SFG(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(_CM, tau, std::string("SFG")) {}
	
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		if(E*(1+z) < 0.6_GeV)
			return pow(1+z, 2-1.5);
		if(E < 0.6_GeV  && 0.6_GeV < E*(1+z))
			return pow(1+z, 2-Gamma_X) * pow(E/0.6_GeV, 1.5 - Gamma_X);
		//if(0.6_GeV <= E)
		return pow(1+z, 2-Gamma_X);
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{
		return pow(E/0.6_GeV, ( E < 0.6_GeV ? -1.5 : -Gamma_X ));
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1;
	}			
	
	
	double LuminosityFunction(const double Luminosity, const double z, const double Gamma) override
	{	// this is equivalent to the rescaled luminosity function, because GammaDistribution = 1
		const double L_solar = 1; // find out
		return infraredLuminosityFunction( exp(alpha* Luminosity/(L_solar*1e10) + beta), z, Gamma)/alpha; // check
	}

protected:
	double infraredLuminosityFunction(const double Luminosity, const double z, const double Gamma) // check
	{
		double L_IR, Phi_s;
		if(z < z_L) 	
		{
			L_IR = L_IR_0 * pow(1+z, k_L1);
			Phi_s = Phi_s_0 * pow(1+z, k_Phi1);
		}	else  
		{
			L_IR = L_IR_0 * pow(1+z_L, k_L1 - k_L2) * pow(1+z, k_L2);
			Phi_s = Phi_s_0 * pow(1+z_Phi, k_Phi1 - k_Phi2) * pow(1+z, k_Phi2);
		}
		return Phi_s * pow(Luminosity/L_IR, 1-alpha) * exp(-pow(log(1.+ Luminosity/L_IR), 2) /(2*sigma*sigma));
	}

};

class NGSFG : public SFG
{
public:
	NGSFG(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(_CM, tau, std::string("NG SFG"))
	{
		alpha = 1.0;
		beta = 39.28;
		sigma = 0.5;
		k_L1 = 4.49;
		k_L2 = 0.;
		z_L = 1.1;
		k_Phi1 = -0.54;
		k_Phi2 = -7.13;
		z_Phi = 0.53;
		Gamma_X = 2.7;
		L_IR_0 = exp(9.78);
		Phi_s_0 = exp(-2.12);
	}
	
};

class SBSFG : public SFG
{
public:
	SBSFG(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(_CM, tau, std::string("SB SFG"))
	{
		alpha = 1.0;
		beta = 39.28;
		sigma = 0.35;
		k_L1 = 1.96;
		k_L2 = 0.;
		z_L = 1.1;
		k_Phi1 = 3.97;
		k_Phi2 = -1.06;
		z_Phi = 1.1;
		Gamma_X = 2.2;
		L_IR_0 = exp(11.17);
		Phi_s_0 = exp(-4.46);
	}
	
};



#endif
