#ifndef SFG_H
#define SFG_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>

class SFG : public AstrophysicalSource
{
protected:
	double gamma_X;
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
	
public:
	SFG(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : AstrophysicalSource(_CM, tau, std::string("SFG"))
	{
		zBounds.first = 1; zBounds.second = 6;
		GammaBounds.first = 0; GammaBounds.second = 0;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		if(E*(1+z) < 0.6_GeV)
			return powf(1+z, 2-1.5);
		if(E < 0.6_GeV  && 0.6_GeV < E*(1+z))
			return powf(1+z, 2-gamma_X) * powf(E/0.6_GeV, 1.5 - gamma_X);
		//if(0.6_GeV <= E)
		return powf(1+z, 2-gamma_X);
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{
		return powf(E/0.6_GeV, ( E < 0.6_GeV ? -1.5 : -gamma_X ));
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
			L_IR = L_IR_0 * powf(1+z, k_L1);
			Phi_s = Phi_s_0 * powf(1+z, k_Phi1);
		}	else  
		{
			L_IR = L_IR_0 * powf(1+z_L, k_L1 - k_L2) * powf(1+z, k_L2);
			Phi_s = Phi_s_0 * powf(1+z_Phi, k_Phi1 - k_Phi2) * powf(1+z, k_Phi2);
		}
		return Phi_s * powf(Luminosity/L_IR, 1-alpha) * exp(-powf(log(1+ Luminosity/L_IR), 2) /(2*sigma*sigma));
	}

};



#endif
