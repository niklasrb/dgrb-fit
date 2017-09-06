#ifndef HALOMODEL_H
#define HALOMODEL_H

/// This header implements a collection of functions that are defined and used by the Halo Model of Cosmology
/// The idea is that you define a bunch of hardcoded functions and the class calculates the functions that are needed for the intensity,
/// 	APS and cross-correlation. This is done through calculating them on the fly or by interpolating, depending on the specific case


#include <memory>
#include <functional>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <string>
#include <gsl/gsl_integration.h>
#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TSpline.h"
#include "Constants.h"
#include "InterpolationWrapper.h"
#include "CosmologyModel.h"
#include "LinearMatterPowerSpectrum.h"


/// This class provides the halo mass function & other important functions which are needed for the Halo Model
class HaloModel
{
	friend class Benchmark; 
	
protected:
	std::shared_ptr<CosmologyModel> CM;					// A pointer to a class modelling a cosmology
	std::shared_ptr<LinearMatterPowerSpectrum> Plin;	// A pointer to a class modelling the linear matter power spectrum
	
	const double CriticalOverdensity = 1.686; 			// Critical Overdensity for spherical collapse
	const double VirialOverdensity = 200.;//18*M_PI*M_PI;		// Or 200?
	Bounds MBounds;										// The mass bounds this Halo Model will modell
	Bounds zBounds;										// redshift bounds this Halo Model will modell
	Bounds kBounds;										// k bounds
	
	double TopHatWindow(const double kR);				// W commonly
	double f_ST(const double v);						// Parametrization of halo mass function by Sheth et al
	double ConcentrationParameter(const double M);		// c_vir or c_200
	double SubhaloBoostFactor(const double M);			// b_sub
	double NFWHaloDensityProfile(const double r, const double M, const double z);	// the Navarro Frenk White profile
public:
	double DMDensity(const double z);					// The dark matter density at a specific redshift
	
	// These routines calculate the std::functions defined below
	void CalculateMatterFluctuationVariance(const std::vector<double>& R, const std::vector<double>& z, std::ostream& out);			
	void CalculateLinearHaloBias(std::vector<double>& M, std::vector<double>& z, std::ostream& out);
	void CalculateHaloMassFunction(const std::vector<double>& M, const std::vector<double>& z, std::ostream& out);
	void CalculateNFWHaloDensityProfileFT(const std::vector<double>& k, const std::vector<double>& M, const std::vector<double>& z, std::ostream& out); //void CalculateNFWHaloDensityProfileFT(const std::vector<double>& kr_Grid, const std::vector<double>& cGrid);
	void CalculateSourceDensitySubhaloBoostFT(std::vector<double>& M, std::vector<double>& k, std::vector<double>& z, std::ostream& out);
	void CalculateClumpingFactor(const std::vector<double>& M, const std::vector<double>& z, std::ostream& out);
	
private:
	struct FTIntegrandParameters;			// Needed for integration
	
public:
	// These functions will be calculated at runtime
	std::function<double(const double M, const double z)> HaloMassFunction = 0;
	std::function<double(const double R, const double z)> MatterFluctuationVariance = 0;
	std::function<double(const double M, const double z)> LinearHaloBias = 0;
	std::function<double(const double k, const double M,  const double z)> NFWHaloDensityProfileFT = 0;
	
	std::function<double(const double k, const double M, const double z)> SourceDensitySubhaloBoostFT = 0;
	std::function<double(const double z)> ClumpingFactor = 0;
	
public:

	HaloModel(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<LinearMatterPowerSpectrum> Plin, Bounds MBounds, Bounds zBounds, Bounds kBounds);
	void Init(unsigned int MLen, unsigned int zLen, unsigned int kLen, std::string file);	 // This guides the calculation of the std::functions
	void Load(std::string file);
	
	void Plot(std::string file);
	

};

HaloModel::HaloModel(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<LinearMatterPowerSpectrum> Plin, Bounds MBounds, Bounds zBounds, Bounds kBounds) : CM(CM), Plin(Plin), MBounds(MBounds) , zBounds(zBounds), kBounds(kBounds)
{
	
}

void HaloModel::Init(unsigned int MLen, unsigned int zLen, unsigned int kLen, std::string file= "")
{
	assert(MLen >=2);	assert(zLen >=2);   assert(kLen >= 2);
	std::fstream f;
	if(!file.empty()) f = std::fstream(file, f.out);
	else f.setstate(std::ios_base::failbit);
	
	std::vector<double> M; M.resize(MLen);
	for(unsigned int i = 0; i < MLen; i++) 
		M.at(i) = exp( log(MBounds.first) + i *(log(MBounds.second) - log(MBounds.first))/(MLen-1.));
	
	std::vector<double> z; z.resize(zLen);
	for(unsigned int i = 0; i < zLen; i++) 
		z.at(i) = exp(log(zBounds.first) + i *(log(zBounds.second) - log(zBounds.first))/(zLen-1.));

	std::vector<double> k; k.resize(kLen);
	for(unsigned int i = 0; i < k.size(); i++)
		k.at(i) = exp(log(kBounds.first) + i *(log(kBounds.second) - log(kBounds.first))/(kLen-1.));
	
	std::vector<double> R; R.resize(M.size());
	for(unsigned int i= 0; i < M.size(); i++)
	{
		R.at(i) = pow(3*M.at(i) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.);
	}
	
	CalculateMatterFluctuationVariance(R, z, f);
	CalculateLinearHaloBias(M, z, f);
	CalculateHaloMassFunction(M, z, f);
	CalculateNFWHaloDensityProfileFT(k, M, z, f);
	CalculateSourceDensitySubhaloBoostFT(M, k, z, f);
	CalculateClumpingFactor(M, z, f);
}



/// Top Hat Window function, commonly W
/// FT of a Top Hat function
/// Implies asymmetric Fourier convention!
double HaloModel::TopHatWindow(const double kR)
{
	return 3.*(sin(kR) - kR*cos(kR))/pow(kR, 3.);
}

/// Calculates the MFV on a grid in R and z and uses the interpolator for the std::func
void HaloModel::CalculateMatterFluctuationVariance(const std::vector<double>& R, const std::vector<double>& z, std::ostream& out)
{
	TF1* Integrand = new TF1("Integrand k^2 * P_lin * W^2  /2pi^2",
											[this] (double* args, double* params) // args[0]: log(k)  params[0]:R  params[1]:z
											{ 	//std::cout << exp(3*args[0]) << '\t' << params[0] << '\t' <<  pow(this->TopHatWindow(exp(args[0])*params[0]) / M_PI, 2) << '\t' << (*Plin)(exp(args[0]), params[1]) << std::endl;
												return exp(3*args[0])* pow(this->TopHatWindow(exp(args[0])*params[0]) / M_PI, 2) * (*Plin)(exp(args[0]), params[1]) /2.; },
											log(this->kBounds.first), log(this->kBounds.second), 2);
	auto MFVSpline = std::make_shared<gsl2DInterpolationWrapper>(R.data(), R.size(), z.data(), z.size());
	
	for(unsigned int i = 0; i < R.size(); i++)
	{
		for(unsigned int j = 0; j < z.size(); j++)
		{
			Integrand->SetParameters(R.at(i), z.at(j));
			MFVSpline->Val(i, j) = Integrand->Integral(log(this->kBounds.first), log(this->kBounds.second), 1e-3);
		}
	}
	
	MFVSpline->Initialize();
	
	std::cout << "MFV: "; MFVSpline->print();
	if(out.good()) MFVSpline->Save(out);
	MatterFluctuationVariance = [MFVSpline] (const double R, const double z)   // more redshift evolution ?
												{	return MFVSpline->Eval(R, z); } ; 
}

/// Part of the Sheth et al Halo Mass Function fit
double HaloModel::f_ST(const double v)
{
	const double A = 0.322;  const double p = 0.3;  const double a = 0.707;
	return A* sqrt(2*a*v/M_PI) * (1. + pow(a*v, -p)) * exp(-a*v/2.);
}

/// Calculates the linear Halo Bias b(m, z)
/// Needs the Matter fluctuation variance
void HaloModel::CalculateLinearHaloBias(std::vector<double>& M, std::vector<double>& z, std::ostream& out)
{
	const double a = 0.707;  const double p = 0.3; double v;
	
	auto LHBSpline = std::make_shared<gsl2DInterpolationWrapper>(M.data(), M.size(), z.data(), z.size());
	
	for(unsigned int i = 0; i < M.size(); i++)
	{
		for(unsigned int j = 0; j < z.size(); j++)
		{
			v = pow(CriticalOverdensity, 2.)/MatterFluctuationVariance(pow(3*M.at(i) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.)*1._Mpc, z.at(j));
			LHBSpline->Val(i, j) = 1. + (a*v - 1.)/CriticalOverdensity + 2.*p/(CriticalOverdensity * (1. + pow(a*v,p)));
		}
	}
	LHBSpline->Initialize();
	
	std::cout << "LHB = "; LHBSpline->print();
	if(out.good()) LHBSpline->Save(out);
	LinearHaloBias = [LHBSpline] (const double M, const double z)
						{  return LHBSpline->Eval(M, z);  };
					
	/*LinearHaloBias = [this, a, p] (const double m, const double z)
					{ const double v = pow(CriticalOverdensity, 2)/MatterFluctuationVariance(pow(3*m / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.), z); 
						//std::cout << "LHB: " << pow(3*m / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.) << '\t' << MatterFluctuationVariance(pow(3*m / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.), z) << '\t' << v << std::endl;
						return 1. + (a*v - 1.)/CriticalOverdensity + 2.*p/(CriticalOverdensity * (1. + pow(a*v,p))); } ;*/
}

/// Calculates Halo Mass Function on a grid and then interpolates using gsl
/// Needs the Matter Fluctuation Variance
void HaloModel::CalculateHaloMassFunction(const std::vector<double>& M, const std::vector<double>& z, std::ostream& out)
{
	auto HMFSpline = std::make_shared<gsl2DInterpolationWrapper>(M.data(), M.size(), z.data(), z.size());
	
	double R =0;			// useful for computation
	std::vector<double> Sigma; Sigma.resize(M.size());
	std::vector<double> logM; logM.resize(M.size());
	std::vector<double> logSigma; logSigma.resize(M.size());
	gsl1DInterpolationWrapper* LogSigmaSpline;				// because we need to interpolate the derivative
	
	for(unsigned int i = 0; i < z.size(); i++)
	{ // Interpolate log(sigma) dependent on log(m) for each z individually and use the derivative
		
		for(unsigned int j = 0; j < M.size(); j++) 
		{
			R = pow(3*M.at(j) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.)*1._Mpc;
			Sigma.at(j) =  sqrt(MatterFluctuationVariance(R, z.at(i))) ;  
			logM.at(j) = log(M.at(j));
			logSigma.at(j) = log(Sigma.at(j));
			//std::cout << "R = " << R << "  sigma= " << Sigma.at(j) << std::endl;
		}
		
		LogSigmaSpline = new gsl1DInterpolationWrapper((double*) logM.data(), logM.size(), (double*) logSigma.data());
		for(unsigned int j = 0; j < M.size(); j++) 
		{
			//std::cout << pow(CriticalOverdensity/Sigma.at(j), 2.) << '\t' << f_ST(pow(CriticalOverdensity/Sigma.at(j), 2.)) << '\t' << std::abs(LogSigmaSpline->Derivative(log(M.at(j)))) << std::endl;
			HMFSpline->Val(j, i) = CM->CriticalDensity*CM->O_m/pow(M.at(j), 2.) * f_ST(pow(CriticalOverdensity/Sigma.at(j), 2.)) * std::abs(LogSigmaSpline->Derivative(logM.at(j)));
		}
		
		delete LogSigmaSpline;
	}
	HMFSpline->Initialize();
	std::cout << "HMF: ";	HMFSpline->print();
	if(out.good()) HMFSpline->Save(out);
	HaloMassFunction = [HMFSpline] (const double M, const double z)
												{ 	return HMFSpline->Eval(M, z); };
}

/// Calculates the concentration parameter c_vir or c_200 for a given mass 
double HaloModel::ConcentrationParameter(const double M)
{
	const double c[] = { 37.5153, -1.5093, 1.636e-2 , 3.66e-4 , -2.89237e-5 , 5.32e-7 };
	double sum = 0;
	for(unsigned int i = 0; i < 6; i++) sum += c[i] * pow(log(M/1._M_solar*h), i);
	return sum;
}

/// Calculates the subhalo boost factor b_sub for a given mass
double HaloModel::SubhaloBoostFactor(const double M)
{
	const double b[] = { -0.442, 0.0796, -0.0025, 4.77e-6 , 4.77e-6 , -9.69e-8};
	double sum = 0;
	for(unsigned int i = 0; i < 6; i++) sum += b[i] * pow(log(M/1._M_solar), i);
	return exp(sum*log(10.));
}

/// Dark Matter density for a givev redshift
double HaloModel::DMDensity(const double z)
{
	return CM->O_dm*CM->CriticalDensity*pow(1.+z, 3.);  // check
}

struct HaloModel::FTIntegrandParameters
{
	double k;  double M;  double z;
	std::function<double(const double r, const double M, const double z)> f;
};

/// Calculates the Fourier Transform of the NFW Density profile
void HaloModel::CalculateNFWHaloDensityProfileFT(const std::vector<double>& k, const std::vector<double>& M, const std::vector<double>& z, std::ostream& out)
{	
	// we need a c grid
	std::vector<double> c; c.resize(2*M.size());
	for(unsigned int i = 0; i < M.size(); i++)
	{
		c.at(2*i) = ConcentrationParameter(M.at(i));
		c.at(2*i+1) = c.at(2*i)/(1.+z[z.size()-1]);
	}
	std::sort(c.begin(), c.end());	
	
	// and a k*r grid is useful
	const double r_vir_min = pow(3*M[0] /(4*M_PI*VirialOverdensity*DMDensity(z[z.size()-1])), 1./3.)*1._Mpc; const double r_s_min = r_vir_min/c[c.size()-1]; 
	const double r_vir_max = pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[0])), 1./3.)*1._Mpc;  const double r_s_max = r_vir_max/c[0];
	Bounds krBounds; krBounds.first = r_s_min*k[0];	krBounds.second = r_s_max*k[k.size()-1];
	std::cout << "krBounds: " << krBounds.first << " - " << krBounds.second << std::endl;
	std::vector<double> kr; kr.resize(2*k.size());
	for(unsigned int i = 0; i < kr.size(); i++)
			kr[i] = exp(log(krBounds.first) + i*(log(krBounds.second) - log(krBounds.first))/(kr.size()-1));
	
	// because of the high oscillations use qawo algorithm
	// https://www.gnu.org/software/gsl/manual/html_node/QAWO-adaptive-integration-for-oscillatory-functions.html
	auto FTIntSpline = std::make_shared<gsl2DInterpolationWrapper>(kr.data(), kr.size(), c.data(), c.size());	
	int ws = 1e4;
	auto workspace = gsl_integration_workspace_alloc(ws);
	auto qawotable = gsl_integration_qawo_table_alloc(0, 0, GSL_INTEG_SINE, std::min(200, ws));
	
	gsl_function Integrand;
	Integrand.function = [](double u, void* params) { return 1./pow(1.+u, 2);  };		// integral is sin(u*kr)/(1+u)^2
	double res = 0, abserr = 0;
	
	for(unsigned int i = 0; i < kr.size(); i++)
	{
		gsl_integration_qawo_table_set(qawotable, kr.at(i), 0, GSL_INTEG_SINE);	// sin(kr * u)
		for(unsigned int j = 0; j < c.size(); j++)
		{
			gsl_integration_qawo_table_set_length(qawotable, c.at(j) - (j > 0 ? c.at(j-1) : 0.));	// integrate in small intervals from c[j-1] to c[j]
			gsl_integration_qawo(&Integrand, (j>0 ? c.at(j-1) : 0.), 0, 1e-4, ws, workspace, qawotable, &res, &abserr);
			FTIntSpline->Val(i, j) = res;
			if(j > 0) FTIntSpline->Val(i, j) += FTIntSpline->Val(i, j-1);
		}
	}	
	gsl_integration_qawo_table_free(qawotable);// clean up	
	gsl_integration_workspace_free(workspace);
	
	FTIntSpline->Initialize();
	std::cout << "FTIntSpline: "; FTIntSpline->print();
	
	auto NFWHDPFTSpline = std::make_shared<Interpolation3DWrapper>(k.data(), k.size(), M.data(), M.size(), z.data(), z.size());	
	
	double c_z, r_vir, r_s ;
	for(unsigned int l = 0; l < z.size(); l++)
	{
		//std::cout << "CalculateNFWHaloDensityProfileFT second loop l = " << l << std::endl;
		for(unsigned int j = 0; j < M.size(); j++)
		{
			c_z = ConcentrationParameter(M.at(j))/(1.+z.at(l));		// don't use c grid because of the std::sort
			r_vir = pow(3.*M.at(j) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))), 1./3.)*1._Mpc;
			r_s = r_vir/c_z;
			for(unsigned int i = 0; i < k.size(); i++)
			{
				//std::cout << i << ", " << j << ", " << l << ": k*r_s = " << k.at(i)*r_s  <<std::endl;
				NFWHDPFTSpline->Val( i, j, l) = M.at(j)*(log(1.+c_z) - c_z/(1.+c_z))/(k.at(i)*r_s*DMDensity(z.at(l))) * FTIntSpline->Eval(k.at(i)*r_s, c_z) ;
			}
		}
	}
	if(out.good()) NFWHDPFTSpline->Save(out);
	NFWHaloDensityProfileFT = [NFWHDPFTSpline] (const double k, const double M, const double z)
												{ /*std::cout << "NFWHDPFTSpline: " << k << ", " << M << ", " << z << std::endl;*/ return NFWHDPFTSpline->Eval(k, M, z); } ; 
	/*
	Bounds rBounds; rBounds.first = pow(3*M[0] /(4*M_PI*VirialOverdensity*DMDensity(z[z.size()-1])), 1./3.);
	rBounds.second = pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[0])), 1./3.);
							
	auto NFWHDPFTSpline = std::make_shared<Interpolation3DWrapper>(k.data(), k.size(), M.data(), M.size(), z.data(), z.size());
	auto cquadws = gsl_integration_cquad_workspace_alloc(100);
	gsl_function Integrand;
	
	Integrand.function = [](double log_r, void* params) //   (FTIntegrandParameters*)params
							{ auto p = (FTIntegrandParameters*)params; return p->f(exp(log_r), p->M, p->z) * exp(2*log_r) * sin(p->k * exp(log_r))/p->k; };
							
	FTIntegrandParameters params;
	Integrand.params = &params;
	
	params.f = [this] (const double r, const double M, const double z) { return NFWHaloDensityProfile(r, M, z)/DMDensity(z); };
	
	double res=0, abserr = 0;
	for(unsigned int i = 0; i < k.size(); i++)
	{
		params.k = k.at(i);
		for(unsigned int j = 0; j < M.size(); j++)
		{
			params.M = M.at(j);
			for(unsigned int l = 0; l < z.size(); l++)
			{
				params.z = z.at(l);
				gsl_integration_cquad(&Integrand, log(1e-30), log(3.*M.at(j) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))))/3., 0, 1e-3, cquadws, &res, &abserr, NULL);
				NFWHDPFTSpline->Val(i, j, l) =  res;
			}
		}
	}
	gsl_integration_cquad_workspace_free(cquadws);
	if(out.good()) NFWHDPFTSpline->Save(out);
	NFWHaloDensityProfileFT = [NFWHDPFTSpline] (const double k, const double M, const double z)
												{ return NFWHDPFTSpline->Eval(k, M, z); } ;  */
}


double HaloModel::NFWHaloDensityProfile(const double r, const double M, const double z)
{
	const double c = ConcentrationParameter(M);
	const double r_vir = pow(3*M /(4*M_PI*VirialOverdensity*DMDensity(z)), 1./3.)*1._Mpc;
	const double r_s = r_vir/c;
	const double rho_s = M/(4*M_PI*pow(r_s,3)) * (log(1.+c) - c/(1.+c));
	return rho_s*r_s/(r*pow(1.+ r/r_s, 2));
}


///Calculates the Fourier transform of the Source Density with a subhalo boost
void HaloModel::CalculateSourceDensitySubhaloBoostFT(std::vector<double>& M, std::vector<double>& k, std::vector<double>& z, std::ostream& out)
{
	// we need a different r grid
	std::vector<double> r; r.resize(M.size());
	Bounds rBounds; rBounds.first = pow(3*M[0] /(4*M_PI*VirialOverdensity*DMDensity(z[z.size()-1])), 1./3.)*1e-1;
	rBounds.second = pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[0])), 1./3.)*1e1;
	
	for(unsigned int i = 0; i < r.size(); i++) r[i] = exp( log(rBounds.first) + i*(log(rBounds.second) - log(rBounds.first))/(r.size()-1.) );
	
	
	auto NFWIntegrand = std::make_shared<TF1>("NFW halo density ^2", 
												[this] ( double* args, double* params) // args[0]: log(r)  params[0]: M  params[1]: z
												{	const double nfw = NFWHaloDensityProfile(exp(args[0]), params[0], params[1]);
													//std::cout << exp(args[0]) << '\t' << nfw << '\t' << DMDensity(params[1]) << '\t' << pow(nfw/DMDensity(params[1]), 2.) << std::endl;
													return exp(args[0])*pow(nfw/DMDensity(params[1]), 2.);	},
													log(rBounds.first),  log(rBounds.second), 2);	
	
	auto SourceDensityWithSubhaloBoost = std::make_shared<Interpolation3DWrapper>(r.data(), r.size(), M.data(), M.size(), z.data(), z.size());
	double NFWsqIntegrated = 0;
	
	for(unsigned int i = 0; i < M.size(); i++)
	{
		//std::cout << "CalculateSourceDensitySubhaloBoostFT first loop i = " << i << std::endl;
		for(unsigned int l = 0; l < z.size(); l++)  
		{
			NFWIntegrand->SetParameters(M.at(i), z.at(l));
			//std::cout << 3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(l))) << '\t' << log(3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(l)))) << std::endl;
			NFWsqIntegrated = NFWIntegrand->Integral(log(rBounds.first), log(3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(l))))/3., 1e-4);
			
			for(unsigned int j = 0; j < r.size(); j++)
			{				
				SourceDensityWithSubhaloBoost->Val(j, i, l) = pow(NFWHaloDensityProfile(r.at(j), M.at(i), z.at(l))/DMDensity(z.at(l)), 2) + 
												SubhaloBoostFactor(M.at(i))* NFWHaloDensityProfile(r.at(j), M.at(i), z.at(l))/M.at(i) * NFWsqIntegrated;
				//std::cout << M.at(i) << '\t' << r.at(j) << '\t' << z.at(l) << '\t' << SourceDensityWithSubhaloBoost->Val(j, i, l) << std::endl;
			}																													// Integrate to r_vir
		}
	}
	
	auto SDSBSpline = std::make_shared<Interpolation3DWrapper>(M.data(), M.size(), k.data(), k.size(), z.data(), z.size());
	
	auto cquadws = gsl_integration_cquad_workspace_alloc(100);
	gsl_function Integrand;
	
	
	Integrand.function = [](double log_r, void* params) // params[0]: k   params[1]: M  params[2]: z     (FTIntegrandParameters*)params
							{ auto p = (FTIntegrandParameters*)params; return p->f(exp(log_r), p->M, p->z) * exp(2*log_r) * sin(p->k * exp(log_r))/p->k; };
							
	FTIntegrandParameters params;
	Integrand.params = &params;
	
	params.f = [SourceDensityWithSubhaloBoost] (const double r, const double M, const double z) { return SourceDensityWithSubhaloBoost->Eval(r, M, z); };
	
	double res=0, abserr = 0;
	for(unsigned int i = 0; i < M.size(); i++)
	{
		//std::cout << "CalculateSourceDensitySubhaloBoostFT second loop i = " << i << std::endl;
		params.M = M.at(i);
		for(unsigned int j = 0; j < k.size(); j++)
		{
			params.k = k.at(j);
			for(unsigned int l = 0; l < z.size(); l++)
			{
				params.z = z.at(l);
				gsl_integration_cquad(&Integrand, log(rBounds.first*1e1), log(3.*M.at(i) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))))/3., 0, 1e-4, cquadws, &res, &abserr, NULL);
				//std::cout << "(" << i << ", " << j << ", " << l << "): " << res << std::endl;
				SDSBSpline->Val(i, j, l) = res;
			}
		}
	}
	gsl_integration_cquad_workspace_free(cquadws);
	
	if(out.good()) SDSBSpline->Save(out);
	SourceDensitySubhaloBoostFT = [SDSBSpline] ( const double k, const double M, const double z)
									{  return SDSBSpline->Eval(M, k, z); };
}

/*
	auto FTIntegrand = std::make_shared<TF1>("Integrand f(r)*r*sin(kr)/k for FT",
										[SourceDensityWithSubhaloBoost] (double* args, double* params) //args[0]:log(r)  params[0]:k  params[1]: M  params[2]: z
										{	return SourceDensityWithSubhaloBoost->Eval(exp(args[0]), params[1], params[2])*exp(2*args[0]) * sin(params[0]*exp(args[0]))/params[0]; },
										log(rBounds.first), log(rBounds.second), 3);
	
	auto SDSBSpline = std::make_shared<Interpolation3DWrapper>(M.data(), M.size(), k.data(), k.size(), z.data(), z.size());
	
	for(unsigned int i = 0; i < M.size(); i++)
	{
		for(unsigned int j = 0; j < k.size(); j++)
		{
			for(unsigned int l = 0; l < z.size(); l++)
			{
				FTIntegrand->SetParameters(k.at(j), M.at(i), z.at(l));
				//std::cout << k.at(j) << '\t' << M.at(i) << '\t' << z.at(l) << '\t' << pow(3.*M.at(i) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))), 1./3.) << std::endl;
				SDSBSpline->Val(i, j, l) = 4*M_PI * FTIntegrand->Integral(log(rBounds.first), log(3.*M.at(i) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))))/3. , 1e-4);
				//std::cout << SDSBSpline->Val(i, j, l) << std::endl;
			}
		}
	}*/

/// Calculates the clumping factor 
void HaloModel::CalculateClumpingFactor(const std::vector<double>& M, const std::vector<double>& z, std::ostream& out)
{
	Bounds rBounds; rBounds.first = pow(3*M[0] /(4*M_PI*VirialOverdensity*DMDensity(z[z.size()-1])), 1./3.)*1e-1;
	rBounds.second = pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[0])), 1./3.)*1e1;
	
	auto NFWIntegrand = std::make_shared<TF1>("NFW halo density ^2", 
												[this] ( double* args, double* params) // args[0]: log(r)  params[0]: M  params[1]: z
												{	const double nfw = NFWHaloDensityProfile(exp(args[0]), params[0], params[1]);
													return exp(args[0])*pow(nfw/DMDensity(params[1]), 2.);	},
													log(rBounds.first),  log(rBounds.second), 2);	
													
	auto NFWsqIntegrated = std::make_shared<gsl2DInterpolationWrapper>(M.data(), M.size(), z.data(), z.size());
	for(unsigned int i = 0; i < M.size(); i++)
	{
		for(unsigned int j = 0; j < z.size(); j++)
		{
			NFWIntegrand->SetParameters(M.at(i), z.at(j));
			NFWsqIntegrated->Val(i, j) = NFWIntegrand->Integral(log(rBounds.first),  log(3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(j))))/3., 1e-4);
		}
	}
	NFWsqIntegrated->Initialize();
	std::cout << "NFWsq "; NFWsqIntegrated->print();
	
	
	auto Integrand = std::make_shared<TF1>("Integrand dn/dm (1+b_sub) Integral(NFW^2)",
											[NFWsqIntegrated, this] (double* args, double* params) // args[0]: log(M)  params[0]: z
											{	//double NFW = NFWsqIntegrated->Eval(exp(args[0]), params[0]);
												//std::cout << exp(args[0]) << '\t' << HaloMassFunction(exp(args[0]), params[0]) << '\t' << (1.+ SubhaloBoostFactor(exp(args[0]))) << '\t' << NFW << std::endl;
												return exp(args[0])*HaloMassFunction(exp(args[0]), params[0])*(1.+ SubhaloBoostFactor(exp(args[0])))*NFWsqIntegrated->Eval(exp(args[0]), params[0]); }, 
												log(MBounds.first), log(MBounds.second), 1);
	std::vector<double> CF;  CF.resize(z.size());
	for(unsigned int i = 0; i < z.size(); i++)
	{
		Integrand->SetParameters(z.at(i), 0);
		CF.at(i) = Integrand->Integral(log(MBounds.first), log(MBounds.second), 1e-4);
	}
	auto CFSpline = std::make_shared<gsl1DInterpolationWrapper>( z.data(), z.size(), CF.data());
	std::cout << "CF: "; CFSpline->print();
	if(out.good()) CFSpline->Save(out);
	ClumpingFactor = [CFSpline] (const double z) { return CFSpline->Eval(z); }; 
	/* ClumpingFactor = [NFWIntegrand, Integrand, this] (const double z) 
												{ Integrand->SetParameters(z, 0);
													return Integrand->Integral(MBounds.first, MBounds.second, 1e-4); } ; */
}

void HaloModel::Load(std::string file)
{
	std::fstream f(file, f.in);
	assert(f.is_open());
	
	auto MFVSpline = std::make_shared<gsl2DInterpolationWrapper>(f);
	MatterFluctuationVariance = [MFVSpline] (const double R, const double z) {	return MFVSpline->Eval(R, z); } ; 
	
	auto LHBSpline = std::make_shared<gsl2DInterpolationWrapper>(f);
	LinearHaloBias = [LHBSpline] (const double M, const double z) { return LHBSpline->Eval(M, z);  };
	
	auto HMFSpline = std::make_shared<gsl2DInterpolationWrapper>(f);
	HaloMassFunction = [HMFSpline] (const double M, const double z)	{ return HMFSpline->Eval(M, z); };
	
	auto NFWHDPFTSpline = std::make_shared<Interpolation3DWrapper>(f);
	NFWHaloDensityProfileFT = [NFWHDPFTSpline] (const double k, const double M, const double z)	{ return NFWHDPFTSpline->Eval(k, M, z); } ;
	
	auto SDSBSpline = std::make_shared<Interpolation3DWrapper>(f);
	SourceDensitySubhaloBoostFT = [SDSBSpline] ( const double k, const double M, const double z) {  return SDSBSpline->Eval(M, k, z); };
									
	auto CFSpline = std::make_shared<gsl1DInterpolationWrapper>(f);
	ClumpingFactor = [CFSpline] (const double z) { return CFSpline->Eval(z); }; 
}

void HaloModel::Plot(std::string file)
{
	auto rootFile = new TFile((file).c_str(), "RECREATE", "HaloModelPlots");
	
	/// Linear Halo Bias
	{
		double* z = new double[4]; z[0]= 0; z[1] = 1; z[2] = 2; z[3] = 3;	
		for(unsigned int i = 0; i < 4; i++)
		{
			TF1 LHB(("LHB2z" + std::to_string(z[i])).c_str(), [this, z, i] (double* args, double* params) { return pow(LinearHaloBias(args[0]*h, z[i]),2.); },
							1e10, 1e14, 0);
			LHB.SetNpx(1e4);
			auto g = new TGraph(&LHB);  g->SetName(("LHB" + std::to_string(z[i])).c_str());
			g->Write();
		}
		delete []z;
	}
	/// Halo Mass Function
	{
		unsigned int n = 30;  double M_min = 1e-6;  double M_max = 1e12;
		double* M = new double[n];
		for(unsigned int i =0; i < n; i++) M[i] = exp( log(M_min) + i*(log(M_max) - log(M_min))/(n-1.));
		double* HMF = new double[n];
		std::vector<int> z = {0, 10, 20, 30};
		for(unsigned int i = 0; i < z.size(); i++)
		{
			for(unsigned int j = 0; j < n; j++) HMF[j] = HaloMassFunction(M[j]*h, z[i]);
			auto g = new TGraph(n, M, HMF); g->SetName(("HMF" + std::to_string(z[i])).c_str());
			g->Write();
		}
		
		delete []M; delete []HMF;
	}
	/// NFWHaloDensityProfileFT
	{
		std::vector<double> M = {1e-5_M_solar, 1._M_solar, 1e10_M_solar};
		std::vector<int> z = {1, 4};
		std::vector<double> k;  k.resize(100); for(unsigned int i = 0; i < k.size(); i++) k[i] = exp(log(kBounds.first) + i*(log(kBounds.second) - log(kBounds.first))/(k.size()-1));
		std::vector<double> NFWHDPFT; NFWHDPFT.resize(k.size());
		for(unsigned int i = 0; i < M.size(); i++)
		{
			for(unsigned int j = 0; j < z.size(); j++)
			{
				//std::cout << i << " " << j << std::endl;
				for(unsigned int l = 0; l < k.size(); l++) NFWHDPFT[l] = NFWHaloDensityProfileFT(k.at(l), M.at(i), z.at(j));
				auto g = new TGraph(k.size(), k.data(), NFWHDPFT.data());
				g->SetName(("NFWHDPFT_" + std::to_string(M.at(i)) + "_" + std::to_string(z.at(j))).c_str());
				g->Write();
			}
		}
	}
	/// SourceDensitySubhaloBoostFT
	{
		std::vector<double> M = {1e-5_M_solar, 1._M_solar, 1e10_M_solar};
		std::vector<int> z = {1, 4};
		std::vector<double> k;  k.resize(100); for(unsigned int i = 0; i < k.size(); i++) k[i] = exp(log(kBounds.first) + i*(log(kBounds.second) - log(kBounds.first))/(k.size()-1));
		std::vector<double> NFWHDPFT; NFWHDPFT.resize(k.size());
		for(unsigned int i = 0; i < M.size(); i++)
		{
			for(unsigned int j = 0; j < z.size(); j++)
			{
				//std::cout << i << " " << j << std::endl;
				for(unsigned int l = 0; l < k.size(); l++) NFWHDPFT[l] = SourceDensitySubhaloBoostFT(k.at(l), M.at(i), z.at(j));
				auto g = new TGraph(k.size(), k.data(), NFWHDPFT.data());
				g->SetName(("SDSBFT_" + std::to_string(M.at(i)) + "_" + std::to_string(z.at(j))).c_str());
				g->Write();
			}
		}
	}
	/// Integrated Halo Mass Function  
	/*Canvas IHMFCanvas("nCanvas", "Integrated Halo Mass Function", 1000, 1000);
	IHMFCanvas().SetLogy(1);
	n = 40;  z = new double[n]; for(unsigned int i = 0; i < n; i++) z[i] = i*(40.)/(n-1.);
	int* Mthres = new int[11]; for(unsigned int i = 0; i < 11; i++) Mthres[i] = 4 + i;
	double* IHMF = new double[n];
	
	TF1 _HMF("HMF", [HM] (double* args, double* params) { return exp(args[0])*HM->HaloMassFunction(h*exp(args[0]), params[0]); }, log(1e4), log(1e18), 1);
	
	for(unsigned int i = 0; i < 11; i++)
	{
		
		for(unsigned int j =0; j < n; j++) 
		{
			_HMF.SetParameters(z[j], 0.);
			IHMF[j] = _HMF.Integral(log(pow(10., Mthres[i])), log(1e18), 1e-4);
		}
		auto IHMFGraph = new TGraph(n, z, IHMF);  IHMFGraph->SetLineColor(i+1);
		IHMFCanvas.AddGraph(IHMFGraph, "$M_{thres}$ = 1e" + std::to_string(Mthres[i]), "L");
	}
	IHMFCanvas.Draw("A");
	IHMFCanvas.SetyLimits(1e-10, 1e10);
	IHMFCanvas.SetxLimits(-1, 40);
	IHMFCanvas().BuildLegend();
	IHMFCanvas().SaveAs("IHMF.jpg");
	
	delete []Mthres; delete []z; */
	
	rootFile->Close();
}

#endif
