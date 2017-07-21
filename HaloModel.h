#ifndef HALOMODEL_H
#define HALOMODEL_H

/// This header implements a collection of functions that are defined and used by the Halo Model of Cosmology
/// The idea is that you define a bunch of hardcoded functions and the class calculates the functions that are needed for the intensity,
/// 	APS and cross-correlation. This is done through calculating them on the fly or by interpolating, depending on the specific case

/// TODO: Calculate clumping factor and FT of source density with subhalo structure
/// 	 to eventually calc the 3D power spectra
#include <memory>
#include <functional>
#include <cmath>
#include <cassert>

#include "TROOT.h"
#include "TF1.h"
#include "TSpline.h"
#include "Constants.h"
#include "gsl2DInterpolationWrapper.h"
#include "CosmologyModel.h"
#include "LinearMatterPowerSpectrum.h"


/// This class provides the halo mass function & other important functions which are needed to calculate the 3D Power Spectrum
class HaloModel
{
	friend class Benchmark;
protected:
	std::shared_ptr<CosmologyModel> CM;					// A pointer to a class modelling a cosmology
	std::shared_ptr<LinearMatterPowerSpectrum> Plin;	// A pointer to a class modelling the linear matter power spectrum
	
	const double CriticalOverdensity = 1.686; 			//Critical Overdensity for spherical collapse
	const double VirialOverdensity = 18*M_PI*M_PI;		// Or 200?
	Bounds MBounds;										// The mass bounds this Halo Model will modell
	Bounds zBounds;										// redshift bounds this Halo Model will modell
	
	double TopHatWindow(const double kR);				// W commonly
	double f_ST(const double v);						// Parametrization of halo mass function by Sheth et al
	double ConcentrationParameter(const double M);		// c_vir or c_200
	double SubhaloBoostFactor(const double M);			// b_sub
	double NFWHaloDensityProfile(const double r, const double M, const double z);	// the Navarro Frenk White profile
	double DMDensity(const double z);					// The dark matter density at a specific redshift
	
	// These routines calculate the std::functions defined below
	void CalculateMatterFluctuationVariance();			
	void CalculateLinearHaloBias();
	void CalculateHaloMassFunction(const std::vector<double>& M, const std::vector<double>& z);
	void CalculateNFWHaloDensityProfileFT();
	void CalculateSourceDensitySubhaloBoostFT();
	void CalculateClumpingFactor(const std::vector<double>& z);
	//void Calculate3DPowerSpectrum();
	
public:
	// These functions will be calculated at runtime
	std::function<double(const double M, const double z)> HaloMassFunction = 0;
	std::function<double(const double R, const double z)> MatterFluctuationVariance = 0;
	std::function<double(const double M, const double z)> LinearHaloBias = 0;
	std::function<double(const double k, const double M,  const double z)> NFWHaloDensityProfileFT = 0;
	
	std::function<double(const double k, const double M, const double z)> SourceDensitySubhaloBoostFT = 0;
	std::function<double(const double z)> ClumpingFactor = 0;
	
public:

	HaloModel(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<LinearMatterPowerSpectrum> Plin, Bounds MBounds, Bounds zBounds);
	void Init(unsigned int MLen, unsigned int zLen); // This guides the calculation of the std::functions
	//void PrepareGalaxyCatalog(GalaxyCatalog* gc);
	

};

HaloModel::HaloModel(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<LinearMatterPowerSpectrum> Plin, Bounds MBounds, Bounds zBounds) : CM(CM), Plin(Plin), MBounds(MBounds) , zBounds(zBounds)
{
	//std::cout << "CM->CriticalDensity: " << CM->CriticalDensity << std::endl;
}

void HaloModel::Init(unsigned int MLen, unsigned int zLen)
{
	assert(MLen >=1);
	assert(zLen >=1);
	std::vector<double> M; M.resize(MLen);
	for(unsigned int i = 0; i < MLen; i++) 
	{
		M.at(i) = exp( log(MBounds.first) + i *(log(MBounds.second) - log(MBounds.first))/MLen);
	}
	std::vector<double> z; z.resize(zLen);
	for(unsigned int i = 0; i < zLen; i++) 
	{
	
		z.at(i) = exp( log(zBounds.first) + i *(log(zBounds.second) - log(zBounds.first))/zLen);
	}
	CalculateMatterFluctuationVariance();
	CalculateLinearHaloBias();
	CalculateHaloMassFunction(M, z);
	CalculateNFWHaloDensityProfileFT();
	CalculateSourceDensitySubhaloBoostFT();
	CalculateClumpingFactor(z);
}

/// Top Hat Window function, commonly W
/// FT of a Top Hat function
/// Implies asymmetric Fourier convention!
double HaloModel::TopHatWindow(const double kR)
{
	//std::cout << kR << '\t' << 2*(sin(kR) - kR*cos(kR)) << '\t' << pow(kR, 3) << '\n';
	return 2*(sin(kR) - kR*cos(kR))/pow(kR, 3);
}

/// Sets the MatterFluctuationVariance Function as one that integrates on every call
void HaloModel::CalculateMatterFluctuationVariance()
{
	auto Integrand = std::make_shared<TF1>("Integrand k^2 * P_lin * W^2  /2pi^2",
											[this] (double* args, double* params) // args[0]: log(k)  params[0]:R  params[1]:z
											{ 	//std::cout << args[0] << '\t' << params[0] << '\t' << pow(args[0]* this->TopHatWindow(args[0]*params[0]) / M_PI, 2) << '\t' << (*Plin)(args[0], params[1]) << std::endl;
												return exp(3*args[0])* pow(this->TopHatWindow(exp(args[0])*params[0]) / M_PI, 2) * (*Plin)(exp(args[0]), params[1]) /2.; },
											0, std::numeric_limits<double>::infinity(), 2);
	MatterFluctuationVariance = [this, Integrand] (const double R, const double z)   // more redshift evolution ?
												{ Integrand->SetParameters(R, z);
													return Integrand->Integral(-3., 4., 1e-4); } ;  // check bounds
}

/// Part of the Sheth et al Halo Mass Function fit
double HaloModel::f_ST(const double v)
{
	const double A = 0.322;  const double p = 0.3;  const double a = 0.707;
	return A* sqrt(2*a*v/M_PI) * (1. + powf(a*v, -p)) * exp(-a*v/2);
}

/// Calculates the linear Halo Bias b(m, z)
/// Needs the Matter fluctuation variance
void HaloModel::CalculateLinearHaloBias()
{
	const double a = 0.707;  const double p = 0.3;
	LinearHaloBias = [this, a, p] (const double m, const double z)
					{ const double v = powf(CM->CriticalDensity, 2)/MatterFluctuationVariance(m, z); 
						return 1+ (a*v - 1)/CM->CriticalDensity + 2*p/(CM->CriticalDensity * (1 + powf(a*v,p))); } ;
}

/// Calculates Halo Mass Function on a grid and then interpolates using gsl
/// Needs the Matter Fluctuation Variance
void HaloModel::CalculateHaloMassFunction(const std::vector<double>& M, const std::vector<double>& z)
{
	double** HMF = new double*[M.size()]; // Holds the HMF values on the grid
	for(unsigned int i= 0; i < M.size(); i++) HMF[i] = new double[z.size()];
	
	double R =0;			// useful for computation
	std::vector<double> Sigma; Sigma.resize(M.size());
	std::vector<double> logM; logM.resize(M.size());
	std::vector<double> logSigma; logSigma.resize(M.size());
	TSpline3* LogSigmaSpline;				// because we need to interpolate the derivative
	
	for(unsigned int i = 0; i < z.size(); i++)
	{ // Interpolate log(sigma) dependent on log(m) for each z individually and use the derivative
		
		for(unsigned int j = 0; j < M.size(); j++) 
		{
			R = powf(3*M.at(j) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.);
			//std::cout << R << std::endl;
			Sigma.at(j) =  sqrt(MatterFluctuationVariance(R, z.at(i))) ;  
			logM.at(j) = log(M.at(j));
			logSigma.at(j) = log(Sigma.at(j));
		}
		
		LogSigmaSpline = new TSpline3("ln(sigma) dependent on log(M) Spline", (double*) logM.data(),(double*) logSigma.data(), M.size());
		for(unsigned int j = 0; j < M.size(); j++) 
		{
			//std::cout << pow(CriticalOverdensity/Sigma.at(j), 2.) << '\t' << f_ST(pow(CriticalOverdensity/Sigma.at(j), 2.)) << '\t' << abs(LogSigmaSpline->Derivative(log(M.at(j)))) << std::endl;
			HMF[j][i] = CriticalOverdensity * CM->O_m * f_ST(pow(CriticalOverdensity/Sigma.at(j), 2.)) * abs(LogSigmaSpline->Derivative(log(M.at(j))));
		}
		
		delete LogSigmaSpline;
	}
	// Now interpolate
	auto HMFSpline = std::make_shared<gsl2DInterpolationWrapper>(M.data(), M.size(), z.data(), z.size(),(const double**) HMF);
	HaloMassFunction = [HMFSpline] (const double M, const double z)
												{ return HMFSpline->Eval(M, z); };
	for(unsigned int i = 0; i < M.size(); i++) delete []HMF[i]; 
	delete []HMF;
}

/// Calculates the concentration parameter c_vir or c_200 for a given mass
double HaloModel::ConcentrationParameter(const double M)
{
	const double c[] = { 37.5153, -1.5093, 1.636e-2 , 3.66e-4 , -2.89237e-5 , 5.32e-7 };
	double sum = 0;
	for(unsigned int i = 0; i < 6; i++) sum += c[i] * powf(log(M/M_solar), i);
	return sum;
}

/// Calculates the subhalo boost factor b_sub for a given mass
double HaloModel::SubhaloBoostFactor(const double M)
{
	const double c[] = { -0.442, 0.0796, -0.0025, 4.77e-6 , 4.77e-6 , -9.69e-8};
	double sum = 0;
	for(unsigned int i = 0; i < 6; i++) sum += c[i] * powf(log(M/M_solar), i);
	return exp(sum);
}

/// rho_DM
double HaloModel::DMDensity(const double z)
{
	return CM->O_dm*CM->CriticalDensity*powf(1+z, 3);  // check
}


/// Calculates the Fourier Transform of the NFW Density profile
void HaloModel::CalculateNFWHaloDensityProfileFT()
{
	auto Integrand = std::make_shared<TF1>("Integrand sin(k*u*r_s)/(1+u^2)", [] (double* args, double *params) // args[0]: u  params[0]: k  params[1]: z  params[2]: r_s
																					{ 	return sin(params[0]* args[0] * params[2]) / (1+args[0]*args[0]); },
																						0,/*range*/ 1e20, 3/*npar*/);			// check range
																						
	NFWHaloDensityProfileFT = [Integrand, this] (const double k, const double M, const double z)
												{ const double c_vir = ConcentrationParameter(M);
													const double r_vir = powf(3*M /(4*M_PI*VirialOverdensity*DMDensity(z)), 1./3.);
													const double r_s = c_vir/r_vir;
													Integrand->SetParameters(k, z, r_s);
													return (log(1+c_vir) - c_vir/(1+c_vir)) * Integrand->Integral(0, c_vir) /(k*r_s); };
}

double HaloModel::NFWHaloDensityProfile(const double r, const double M, const double z)
{
	const double c = ConcentrationParameter(M);
	const double r_vir = pow(3*M /(4*M_PI*VirialOverdensity*DMDensity(z)), 1./3.);
	const double r_s = c/r_vir;
	const double rho_s = M/(4*M_PI*pow(r_s,3)) * (log(1+c) - c/(1+c));
	return rho_s*r_s/(r*pow(1+ r/r_s, 2));
}

///Calculates the Fourier transform of the Source Density with a subhalo boost
void HaloModel::CalculateSourceDensitySubhaloBoostFT()
{
	auto NFWIntegrand = std::make_shared<TF1>("NFW halo density ^2", 
												[this] ( double* args, double* params) // args[0]: r  params[0]: M  params[1]: z
												{
													return 4*M_PI*powf(NFWHaloDensityProfile(args[0], params[0], params[1])/DMDensity(params[1]), 2);	},
													0,  1e99, 2);		// check range
	
	std::function<double(const double, const double, const double)> SourceDensityWithSubhaloBoost = 
							[this, NFWIntegrand] (const double r, const double M, const double z) //args[0]: r params[0]: M  params[1]: z 
							{	NFWIntegrand->SetParameters(M, z);
								return powf(NFWHaloDensityProfile(r, M, z)/DMDensity(z), 2) + 
												SubhaloBoostFactor(M)* NFWHaloDensityProfile(r, M, z)/M * NFWIntegrand->Integral(0, powf(3*M /(4*M_PI*VirialOverdensity*DMDensity(z)), 1./3.), 1e-4); // Integrate to r_vir
							}; 
	
	// TODO: calculate this on a grid to avoid the thousand nested integrals
	
	auto FTIntegrand = std::make_shared<TF1>("Integrand f(r)*r*sin(kr)/k for FT",
										[SourceDensityWithSubhaloBoost] (double* args, double* params) //args[0]:r  params[0]:k  params[1]: M  params[2]: z
										{	return SourceDensityWithSubhaloBoost(args[0], params[1], params[2])*args[0] * sin(params[0]*args[0])/params[0]; },
										0, 1e99, 3);
	
	SourceDensitySubhaloBoostFT = [this, FTIntegrand] ( const double k, const double M, const double z)
									{  FTIntegrand->SetParameters(k, M, z);
										return 4*M_PI * FTIntegrand->Integral(0, powf(3*M /(4*M_PI*VirialOverdensity*DMDensity(z)), 1./3.), 1e-4); };
}

/// Calculates the clumping factor 
void HaloModel::CalculateClumpingFactor(const std::vector<double>& z)
{
	auto NFWIntegrand = std::make_shared<TF1>("NFW halo density ^2", 
												[this] ( double* args, double* params) // args[0]: r  params[0]: M  params[1]: z
												{ 	double ret = 4*M_PI*pow(NFWHaloDensityProfile(args[0], params[0], params[1])/DMDensity(params[1]), 2);
													//std::cout << ret << std::endl;
													return 	ret;},
													0,  1e99, 2);		// check range
													
	auto Integrand = std::make_shared<TF1>("Integrand dn/dm (1+b_sub) Integral(NFW^2)",
											[NFWIntegrand, this] (double* args, double* params) // args[0]: log(M)  params[0]: z
											{	NFWIntegrand->SetParameters(exp(args[0]), params[0]);
												double NFW = NFWIntegrand->Integral(0,  pow(3*exp(args[0]) /(4*M_PI*VirialOverdensity*DMDensity(params[0])), 1./3.), 1e-4);
												std::cout << exp(args[0]) << '\t' << HaloMassFunction(exp(args[0]), params[0]) << '\t' << (1.+ SubhaloBoostFactor(exp(args[0]))) << '\t' << NFW << std::endl;
												return exp(args[0])*HaloMassFunction(exp(args[0]), params[0])*(1.+ SubhaloBoostFactor(exp(args[0])))*std::max(0., NFW); }, // Integrate to r_vir
												log(MBounds.first), log(MBounds.second), 1);
	std::vector<double> CF;  CF.resize(z.size());
	for(unsigned int i = 0; i < z.size(); i++)
	{
		Integrand->SetParameters(z.at(i), 0);
		CF.at(i) = Integrand->Integral(log(MBounds.first), log(MBounds.second), 1e-4);
	}
	auto CFSpline = std::make_shared<gsl1DInterpolationWrapper>( z.data(), z.size(), CF.data());
	CFSpline->print();
	ClumpingFactor = [CFSpline] (const double z) { return CFSpline->Eval(z); }; 
	/* ClumpingFactor = [NFWIntegrand, Integrand, this] (const double z) 
												{ Integrand->SetParameters(z, 0);
													return Integrand->Integral(MBounds.first, MBounds.second, 1e-4); } ; */
}

/// Calculates source density field, its fourier transform and effective galaxy density
/*void HaloModel::PrepareGalaxyCatalog(GalaxyCatalog* gc)
{
	auto EffectiveGalaxyDensityIntegrand = std::make_shared<TF1>((std::string("Integrand dn/dm N_g for ") + gc->name).c_str(),
											[gc, this] (double* args, double* params) // args[0]: M  params[0]:z
											{ return HaloMassFunction(args[0], params[0])*(gc->N_central(args[0]) + gc->N_satellite(args[0])); },
											MBounds.first, MBounds.second, 1);
	gc->EffectiveGalaxyDensity = [EffectiveGalaxyDensityIntegrand, this] (const double z)
											{ return EffectiveGalaxyDensityIntegrand->Integral(MBounds.first, MBounds.second, 1e-4); };
	
	gc->SourceDensity = [gc, this] (const double r, const double M, const double z)
						{ return gc->N_central(M) + gc->N_satellite(M) * NFWHaloDensityProfile(r, M, z) / M; };
	
	gc->SourceDensityFT = [gc, this] (const double k, const double M, const double z)
						{ return gc->N_central(M) + gc->N_satellite(M) * NFWHaloDensityProfileFT(k, M, z) / M; };
}*/

#endif
