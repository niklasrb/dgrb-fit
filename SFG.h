#ifndef SFG_H
#define SFG_H

#include "AstrophysicalSource.h"
#include "Constants.h"

#include <cmath>
#include <algorithm>
#include <tuple>
#include <string>


class SFG;
class SBSFG;
class SFGAGN;
class NGSFG;

/// Star forming Galaxies
class SFG : public AstrophysicalSource
{
protected:
	double Gamma_X;
	double alpha;
	const double beta = 39.28;
	const double alphatransform = 1.17;
	double sigma;
	double k_L1;
	double k_L2;
	double z_L;
	double k_Phi1;
	double k_Phi2;
	double z_Phi;
	double L_IR_0; 
	double Phi_s_0;
	double E_cut = 0.6_GeV;
	double k;
	
	SFG(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau, std::string name) : AstrophysicalSource(CM, tau, name)
	{
		zBounds.first = 0; zBounds.second = 4;
		GammaBounds.first = 0; GammaBounds.second = 0;
		LumBounds.first = 1e20_ergpers; LumBounds.second = 1e52_ergpers;
	}
	
public:
	SFG(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(CM, tau, std::string("SFG")) {}
	
	SFG(std::shared_ptr<SBSFG> SB, std::shared_ptr<NGSFG> NG, std::shared_ptr<SFGAGN> AGN) ;
	
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{
		if(E*(1.+z) < E_cut)
			return pow(1.+z, 2.-1.5);
		if(E < E_cut  && E_cut <= E*(1.+z))
			return pow(1.+z, 2.-Gamma_X) * pow(E/E_cut, 1.5 - Gamma_X);
		return pow(1.+z, 2.-Gamma_X);
	}
	
	double literal_F(const double E, const double z, const double Gamma) override
	{
		return pow(E/E_cut, ( E < E_cut ? -1.5 : -Gamma_X ));
	}
	
	double GammaDistribution(const double Gamma) override
	{
		return 1;
	}			
	
	double LuminosityFunction(const double L, const double z) override
	{
		const double L_IR = 1e10*L_solar  *pow(L/1._ergpers, 1./alphatransform)*pow(10., -beta/alphatransform);
		return k*infraredLuminosityFunction( L_IR, z)/alphatransform;
	}

protected:
	double infraredLuminosityFunction(const double Luminosity, const double z) 
	{
		double L_IR, Phi_s;
		if(z < z_L) 	
			L_IR = L_IR_0 * pow(1.+z, k_L1);
		else  
			L_IR = L_IR_0 * pow(1.+z_L, k_L1 - k_L2) * pow(1+z, k_L2);
		
		if(z < z_Phi)
			Phi_s = Phi_s_0 * pow(1.+z, k_Phi1);
		else
			Phi_s = Phi_s_0 * pow(1.+z_Phi, k_Phi1 - k_Phi2) * pow(1.+z, k_Phi2);
		return Phi_s * pow(Luminosity/L_IR, 1.-alpha) * exp(-pow(log(1.+ Luminosity/L_IR), 2) /(2*sigma*sigma));
	}

};

class NGSFG : public SFG	// normal galaxies that are star forming galaxies
{
public:
	NGSFG(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(CM, tau, std::string("SFG NG"))
	{
		alpha = 1.0;
		sigma = 0.5;
		k_L1 = 4.49;
		k_L2 = 0.;
		z_L = 1.1;
		k_Phi1 = -0.54;
		k_Phi2 = -7.13;
		z_Phi = 0.53;
		Gamma_X = 2.7;
		L_IR_0 = pow(10, 9.78)*L_solar;		
		Phi_s_0 = pow(10, -2.12)/pow(1._Mpc, 3);
		
		k = 1e5;
	}
	
};

class SBSFG : public SFG	// star burst star forming galaxies
{
public:
	SBSFG(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(CM, tau, std::string("SFG SB"))
	{
		alpha = 1.0;
		sigma = 0.35;
		k_L1 = 1.96;
		k_L2 = 0.;
		z_L = 100;
		k_Phi1 = 3.97;
		k_Phi2 = -1.06;
		z_Phi = 1.1;
		Gamma_X = 2.2;
		L_IR_0 = pow(10., 11.17)*L_solar;	
		Phi_s_0 = pow(10., -4.46)/pow(1._Mpc, 3);
		
		
		LumBounds.first = 1e35;
		k = 1e5;
	}
	
};

class SFGAGN : public SFG	// Star forming galaxies with active galactic nuclei
{
private:
	std::shared_ptr<SBSFG> sbsfg;
	std::shared_ptr<NGSFG> ngsfg;
	
	std::vector<std::pair<Bounds, double> > redshift_dependence;
public:
	SFGAGN(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<EBLAbsorbtionCoefficient> tau) : SFG(CM, tau, std::string("SFG AGN"))
	{
		alpha = 1.20;
		sigma = 0.40;
		k_L1 = 3.17;
		k_L2 = 0.;
		z_L = 100;
		k_Phi1 = 0.67;
		k_Phi2 = -3.17;
		z_Phi = 1.1;
		L_IR_0 = pow(10, 10.80)*L_solar;	
		Phi_s_0 = pow(10, -3.2)/pow(1._Mpc, 3);
		
		k = 1e5;
		
		sbsfg = std::make_shared<SBSFG>(CM, tau);
		ngsfg = std::make_shared<NGSFG>(CM, tau);
		redshift_dependence = {std::make_pair(Bounds(0., 0.3), 0.85), std::make_pair(Bounds(0.3, 0.45), 0.91), std::make_pair(Bounds(0.45, 0.6), 0.99),
								std::make_pair(Bounds(0.6, 0.8), 0.87), std::make_pair(Bounds(0.8, 1.0), 0.73), std::make_pair(Bounds(1.0, 1.2), 0.32),
								std::make_pair(Bounds(1.2, 1.7), 0.75), std::make_pair(Bounds(1.7, 2.0), 0.75), std::make_pair(Bounds(2.0, 2.5), 0.19),
								std::make_pair(Bounds(2.5, 3.0), 0.24), std::make_pair(Bounds(3.0, 4.2), 0.28)};
	}
	

	double literal_F(const double E, const double z, const double Gamma) override
	{
		const double sb = sbsfg->literal_F(E, z, Gamma);
		const double ng = ngsfg->literal_F(E, z, Gamma);
		for(unsigned int i = 0; i < redshift_dependence.size(); i++)
		{
			if(redshift_dependence[i].first.first < z && z <= redshift_dependence[i].first.second)  return ng*redshift_dependence[i].second + sb*(1.-redshift_dependence[i].second);
		}
		return 0;
	}
	
	double EnergySpectrumOverK(const double E, const double z, const double Gamma) override
	{
		const double sb = sbsfg->EnergySpectrumOverK(E, z, Gamma);
		const double ng = ngsfg->EnergySpectrumOverK(E, z, Gamma);
		for(unsigned int i = 0; i < redshift_dependence.size(); i++)
		{
			if(redshift_dependence[i].first.first < z && z <= redshift_dependence[i].first.second)  return ng*redshift_dependence[i].second + sb*(1.-redshift_dependence[i].second);
		}
		return 0;
	}
};

/// This constructor cobines all the different SFG populations into one 
/// By adding the results for intensity and APS coeffiecients together
SFG::SFG(std::shared_ptr<SBSFG> SB, std::shared_ptr<NGSFG> NG, std::shared_ptr<SFGAGN> AGN) : SFG(SB->CM, SB->tau, std::string("SFG"))	
{
	assert(SB->Intensity.size() == NG->Intensity.size() && NG->Intensity.size() == AGN->Intensity.size());
	assert(SB->APS->Bin1Size() == NG->APS->Bin1Size() && NG->APS->Bin1Size() == AGN->APS->Bin1Size());
	this->Intensity.resize(SB->Intensity.size());
	for(unsigned int i = 0; i < Intensity.size(); i++) Intensity[i] = SB->Intensity[i] + NG->Intensity[i] + AGN->Intensity[i];
	APS = std::make_shared<AstrophysicalSourceAPS<std::shared_ptr<gsl1DInterpolationWrapper> > >(SB->APS->Bin1Size(), SB->APS->Bin2Size());
	for(unsigned int i = 0; i < SB->APS->Bin1Size(); i++)
	{
		for(unsigned int j = 0; j <= i; j++)
		{
			APS->at(i, j) = std::make_shared<gsl1DInterpolationWrapper>(*(SB->APS->at(i,j)) + *(NG->APS->at(i,j)) + *(AGN->APS->at(i,j)));
			if(i!=j) APS->at(j, i) = APS->at(i, j);
		}
	}
}



#endif
