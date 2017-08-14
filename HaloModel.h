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
#include <algorithm>

#include "TROOT.h"
#include "TF1.h"
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
	double DMDensity(const double z);					// The dark matter density at a specific redshift
	
	// These routines calculate the std::functions defined below
	void CalculateMatterFluctuationVariance(const std::vector<double>& R, const std::vector<double>& z);			
	void CalculateLinearHaloBias(std::vector<double>& M, std::vector<double>& z);
	void CalculateHaloMassFunction(const std::vector<double>& M, const std::vector<double>& z);
	void CalculateNFWHaloDensityProfileFT(const std::vector<double>& kr_Grid, const std::vector<double>& cGrid);
	void CalculateSourceDensitySubhaloBoostFT(std::vector<double>& M, std::vector<double>& k, std::vector<double>& z);
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

	HaloModel(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<LinearMatterPowerSpectrum> Plin, Bounds MBounds, Bounds zBounds, Bounds kBounds);
	void Init(unsigned int MLen, unsigned int zLen, unsigned int kLen);	 // This guides the calculation of the std::functions
	//void PrepareGalaxyCatalog(GalaxyCatalog* gc);
	

};

HaloModel::HaloModel(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<LinearMatterPowerSpectrum> Plin, Bounds MBounds, Bounds zBounds, Bounds kBounds) : CM(CM), Plin(Plin), MBounds(MBounds) , zBounds(zBounds), kBounds(kBounds)
{
	//std::cout << "CM->CriticalDensity: " << CM->CriticalDensity << std::endl;
}

void HaloModel::Init(unsigned int MLen, unsigned int zLen, unsigned int kLen)
{
	assert(MLen >=1);	assert(zLen >=1);   assert(kLen >= 1);
	
	std::vector<double> M; M.resize(MLen);
	for(unsigned int i = 0; i < MLen; i++) 
		M.at(i) = exp( log(MBounds.first) + i *(log(MBounds.second) - log(MBounds.first))/(MLen-1.));
	
	std::vector<double> z; z.resize(zLen);
	for(unsigned int i = 0; i < zLen; i++) 
		z.at(i) = exp(log(zBounds.first) + i *(log(zBounds.second) - log(zBounds.first))/(zLen-1.));

	std::vector<double> k; k.resize(kLen);
	for(unsigned int i = 0; i < k.size(); i++)
		k.at(i) = exp(log(kBounds.first) + i *(log(kBounds.second) - log(kBounds.first))/(kLen-1.));
	
	std::vector<double> c; c.resize(M.size());
	for(unsigned int i = 0; i < c.size(); i++)
		c.at(i) = ConcentrationParameter(M.at(i));
	std::sort(c.begin(), c.end());
	
	std::vector<double> R; R.resize(M.size()/*z.size()*/);
	for(unsigned int i= 0; i < M.size(); i++)
	{
		R.at(i) = pow(3*M.at(i) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.);
	}
	//std::cout << "R[0] = " << R[0] << "   R[" << R.size()-1 << "] = " << R[R.size() - 1] << std::endl;
	
	CalculateMatterFluctuationVariance(R, z);
	CalculateLinearHaloBias(M, z);
	CalculateHaloMassFunction(M, z);
	CalculateNFWHaloDensityProfileFT(k, c);
	CalculateSourceDensitySubhaloBoostFT(M, k, z);
	CalculateClumpingFactor(z);
}

/// Top Hat Window function, commonly W
/// FT of a Top Hat function
/// Implies asymmetric Fourier convention!
double HaloModel::TopHatWindow(const double kR)
{
	//std::cout << kR << '\t' << 2*(sin(kR) - kR*cos(kR)) << '\t' << pow(kR, 3) << '\n';
	return 2.*(sin(kR) - kR*cos(kR))/pow(kR, 3);
}

/// Sets the MatterFluctuationVariance Function as one that integrates on every call
void HaloModel::CalculateMatterFluctuationVariance(const std::vector<double>& R, const std::vector<double>& z)
{
	TF1* Integrand = new TF1("Integrand k^2 * P_lin * W^2  /2pi^2",
											[this] (double* args, double* params) // args[0]: log(k)  params[0]:R  params[1]:z
											{ 	//std::cout << exp(3*args[0]) << '\t' << params[0] << '\t' <<  pow(this->TopHatWindow(exp(args[0])*params[0]) / M_PI, 2) << '\t' << (*Plin)(exp(args[0]), params[1]) << std::endl;
												return exp(3*args[0])* pow(this->TopHatWindow(exp(args[0])*params[0]) / M_PI, 2) * (*Plin)(exp(args[0]), params[1]) /2.; },
											log(this->kBounds.first), log(this->kBounds.second), 2);
	double** MFV = new double*[R.size()];
	for(unsigned int i = 0; i < R.size(); i++) MFV[i] = new double[z.size()];
	for(unsigned int i = 0; i < R.size(); i++)
	{
		for(unsigned int j = 0; j < z.size(); j++)
		{
			Integrand->SetParameters(R.at(i), z.at(j));
			MFV[i][j] = Integrand->Integral(log(this->kBounds.first), log(this->kBounds.second), 1e-3);
		}
	}
	
	auto MFVSpline = std::make_shared<gsl2DInterpolationWrapper>(R.data(), R.size(), z.data(), z.size(), (const double**) MFV);
	
	std::cout << "MFV: "; MFVSpline->print();
	
	MatterFluctuationVariance = [MFVSpline] (const double R, const double z)   // more redshift evolution ?
												{ if(MFVSpline->Eval(R, z) <= 0) std::cout << "MFV: R = " << R << "  z = " << z << " MFV = " << MFVSpline->Eval(R, z) << std::endl;
													return MFVSpline->Eval(R, z); } ; 
}

/// Part of the Sheth et al Halo Mass Function fit
double HaloModel::f_ST(const double v)
{
	const double A = 0.322;  const double p = 0.3;  const double a = 0.707;
	return A* sqrt(2*a*v/M_PI) * (1. + pow(a*v, -p)) * exp(-a*v/2.);
}

/// Calculates the linear Halo Bias b(m, z)
/// Needs the Matter fluctuation variance
void HaloModel::CalculateLinearHaloBias(std::vector<double>& M, std::vector<double>& z)
{
	const double a = 0.707;  const double p = 0.3; double v;
	
	double** LHB = new double*[M.size()];
	for(unsigned int i = 0; i < M.size(); i++) LHB[i] = new double[z.size()];
	
	for(unsigned int i = 0; i < M.size(); i++)
	{
		for(unsigned int j = 0; j < z.size(); j++)
		{
			v = pow(CriticalOverdensity, 2)/MatterFluctuationVariance(pow(3*M.at(i) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.), z.at(j));
			LHB[i][j] = 1. + (a*v - 1.)/CriticalOverdensity + 2.*p/(CriticalOverdensity * (1. + pow(a*v,p)));
		}
	}
	
	auto LHBSpline = std::make_shared<gsl2DInterpolationWrapper>(M.data(), M.size(), z.data(), z.size(), (const double**) LHB);
	std::cout << "LHB = "; LHBSpline->print();
	LinearHaloBias = [LHBSpline] (const double M, const double z)
						{  return LHBSpline->Eval(M, z);  };
					
	for(unsigned int i = 0; i < M.size(); i++) delete [] LHB[i];
	delete []LHB;
	/*LinearHaloBias = [this, a, p] (const double m, const double z)
					{ const double v = pow(CriticalOverdensity, 2)/MatterFluctuationVariance(pow(3*m / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.), z); 
						//std::cout << "LHB: " << pow(3*m / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.) << '\t' << MatterFluctuationVariance(pow(3*m / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.), z) << '\t' << v << std::endl;
						return 1. + (a*v - 1.)/CriticalOverdensity + 2.*p/(CriticalOverdensity * (1. + pow(a*v,p))); } ;*/
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
	gsl1DInterpolationWrapper* LogSigmaSpline;				// because we need to interpolate the derivative
	
	for(unsigned int i = 0; i < z.size(); i++)
	{ // Interpolate log(sigma) dependent on log(m) for each z individually and use the derivative
		
		for(unsigned int j = 0; j < M.size(); j++) 
		{
			R = pow(3*M.at(j) / (4.*M_PI * CM->CriticalDensity * CM->O_m), 1./3.);
			Sigma.at(j) =  sqrt(MatterFluctuationVariance(R, z.at(i))) ;  
			logM.at(j) = log(M.at(j));
			logSigma.at(j) = log(Sigma.at(j));
			//std::cout << "R = " << R << "  sigma= " << Sigma.at(j) << std::endl;
		}
		
		LogSigmaSpline = new gsl1DInterpolationWrapper((double*) logM.data(), logM.size(), (double*) logSigma.data());
		for(unsigned int j = 0; j < M.size(); j++) 
		{
			//std::cout << pow(CriticalOverdensity/Sigma.at(j), 2.) << '\t' << f_ST(pow(CriticalOverdensity/Sigma.at(j), 2.)) << '\t' << std::abs(LogSigmaSpline->Derivative(log(M.at(j)))) << std::endl;
			HMF[j][i] = CM->CriticalDensity*CM->O_m/M.at(j) * f_ST(pow(CriticalOverdensity/Sigma.at(j), 2.)) * std::abs(LogSigmaSpline->Derivative(logM.at(j)));
		}
		
		delete LogSigmaSpline;
	}
	// Now interpolate
	auto HMFSpline = std::make_shared<gsl2DInterpolationWrapper>(M.data(), M.size(), z.data(), z.size(),(const double**) HMF);
	std::cout << "HMF: ";	HMFSpline->print();
	
	HaloMassFunction = [HMFSpline] (const double M, const double z)
												{ //std::cout << M << '\t' << z << '\t' << HMFSpline->Eval(M, z) << std::endl;
													return HMFSpline->Eval(M, z); };
	for(unsigned int i = 0; i < M.size(); i++) delete []HMF[i]; 
	delete []HMF;
}

/// Calculates the concentration parameter c_vir or c_200 for a given mass in M_solar
double HaloModel::ConcentrationParameter(const double M)
{
	const double c[] = { 37.5153, -1.5093, 1.636e-2 , 3.66e-4 , -2.89237e-5 , 5.32e-7 };
	double sum = 0;
	for(unsigned int i = 0; i < 6; i++) sum += c[i] * pow(log(M*h), i);
	return sum;
}

/// Calculates the subhalo boost factor b_sub for a given mass in M_solar
double HaloModel::SubhaloBoostFactor(const double M)
{
	const double b[] = { -0.442, 0.0796, -0.0025, 4.77e-6 , 4.77e-6 , -9.69e-8};
	double sum = 0;
	for(unsigned int i = 0; i < 6; i++) sum += b[i] * pow(log(M), i);
	return exp(sum*log(10.));
}

/// Dark Matter density for a givev redshift
double HaloModel::DMDensity(const double z)
{
	return CM->O_dm*CM->CriticalDensity*pow(1.+z, 3.);  // check
}


/// Calculates the Fourier Transform of the NFW Density profile
void HaloModel::CalculateNFWHaloDensityProfileFT(const std::vector<double>& kr_Grid, const std::vector<double>& cGrid)
{
	auto Integrand = std::make_shared<TF1>("Integrand sin(k*u*r_s)/(1+u^2)", [] (double* args, double *params) // args[0]: u  params[0]: k*r_s
																					{ 	return std::sin(params[0]*args[0]) / (1.+pow(args[0],2)); },
																						0,/*range*/ 1e20, 1/*npar*/);			// check range
	
	double** FTInt = new double*[kr_Grid.size()];
	for(unsigned int i = 0; i < kr_Grid.size(); i++) FTInt[i] = new double[cGrid.size()];
	
	for(unsigned int i = 0; i < kr_Grid.size(); i++)
	{
		Integrand->SetParameters(kr_Grid[i], 0.);
		
		for(unsigned int j = 0; j < cGrid.size(); j++)
		{
			FTInt[i][j] = Integrand->Integral((j > 0 ? cGrid[j-1] : 0), cGrid[j], 1e-4);
			if(j > 0) FTInt[i][j] += FTInt[i][j-1];
		}
	}
	auto FTIntSpline = std::make_shared<gsl2DInterpolationWrapper>(kr_Grid.data(), kr_Grid.size(), cGrid.data(), cGrid.size(), (const double**)FTInt);
	//std::cout << "FTIntSpline: "; FTIntSpline->print();
																						
	NFWHaloDensityProfileFT = [FTIntSpline, this] (const double k, const double M, const double z)
												{ const double c_vir = ConcentrationParameter(M)/(1.+z);
													const double r_vir = pow(3.*M /(4.*M_PI*VirialOverdensity*DMDensity(z)), 1./3.);
													const double r_s = c_vir/r_vir;
													//std::cout << "NFW: " << k << '\t' << M << '\t' << z << '\t' << c_vir << '\t' << r_vir << std::endl;
													return (log(1.+c_vir) - c_vir/(1.+c_vir)) * FTIntSpline->Eval(k*r_s, c_vir) /(k*r_s); };
}

double HaloModel::NFWHaloDensityProfile(const double r, const double M, const double z)
{
	const double c = ConcentrationParameter(M);
	const double r_vir = pow(3*M /(4*M_PI*VirialOverdensity*DMDensity(z)), 1./3.);
	const double r_s = c/r_vir;
	const double rho_s = M/(4*M_PI*pow(r_s,3)) * (log(1.+c) - c/(1.+c));
	//std::cout << "NFW: " << r << '\t' << M << '\t' << c << '\t' << r_vir << '\t' << r_s << '\t' << rho_s << '\t' << rho_s*r_s/(r*pow(1.+ r/r_s, 2)) << std::endl;
	return rho_s*r_s/(r*pow(1.+ r/r_s, 2));
}

///Calculates the Fourier transform of the Source Density with a subhalo boost
void HaloModel::CalculateSourceDensitySubhaloBoostFT(std::vector<double>& M, std::vector<double>& k, std::vector<double>& z)
{
	// we need a different r grid
	std::vector<double> r; r.resize(M.size());
	Bounds rBounds; rBounds.first = pow(3*M[0] /(4*M_PI*VirialOverdensity*DMDensity(z[z.size()-1])), 1./3.)*1e-10;
	rBounds.second = pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[0])), 1./3.);
	//std::cout << pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[z.size()-1])), 1./3.) << '\t' << pow(3*M[M.size()-1] /(4*M_PI*VirialOverdensity*DMDensity(z[0])), 1./3.) <<std::endl;
	
	for(unsigned int i = 0; i < r.size(); i++) r[i] = exp( log(rBounds.first) + i*(log(rBounds.second) - log(rBounds.first))/(r.size()-1.) );
	
	
	auto NFWIntegrand = std::make_shared<TF1>("NFW halo density ^2", 
												[this] ( double* args, double* params) // args[0]: log(r)  params[0]: M  params[1]: z
												{	const double nfw = NFWHaloDensityProfile(exp(args[0]), params[0], params[1]);
													//std::cout << exp(args[0]) << '\t' << nfw << '\t' << DMDensity(params[1]) << '\t' << pow(nfw/DMDensity(params[1]), 2.) << std::endl;
													return exp(args[0])*pow(nfw/DMDensity(params[1]), 2.);	},
													log(rBounds.first),  log(rBounds.second), 2);		// check range
	
	auto SourceDensityWithSubhaloBoost = std::make_shared<Interpolation3DWrapper>(r.data(), r.size(), M.data(), M.size(), z.data(), z.size());
	double NFWsqIntegrated = 0;
	
	for(unsigned int i = 0; i < M.size(); i++)
	{
		for(unsigned int l = 0; l < z.size(); l++)  
		{
			NFWIntegrand->SetParameters(M.at(i), z.at(l));
			std::cout << 3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(l))) << '\t' << log(3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(l)))) << std::endl;
			NFWsqIntegrated = NFWIntegrand->Integral(log(rBounds.first), log(3*M.at(i) /(4*M_PI*VirialOverdensity*DMDensity(z.at(l))))/3., 1e-4);
			
			for(unsigned int j = 0; j < r.size(); j++)
			{				
				SourceDensityWithSubhaloBoost->Val(j, i, l) = pow(NFWHaloDensityProfile(r.at(j), M.at(i), z.at(l))/DMDensity(z.at(l)), 2) + 
												SubhaloBoostFactor(M.at(i))* NFWHaloDensityProfile(r.at(j), M.at(i), z.at(l))/M.at(i) * NFWsqIntegrated;
				std::cout << M.at(i) << '\t' << r.at(j) << '\t' << z.at(l) << '\t' << SourceDensityWithSubhaloBoost->Val(j, i, l) << std::endl;
			}																													// Integrate to r_vir
		}
	}
	//SourceDensityWithSubhaloBoost->print();'
	// TODO: calculate this on a grid to avoid the thousand nested integrals
	
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
				std::cout << k.at(j) << '\t' << M.at(i) << '\t' << z.at(l) << '\t' << pow(3.*M.at(i) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))), 1./3.) << std::endl;
				SDSBSpline->Val(i, j, l) = 4*M_PI * FTIntegrand->Integral(log(rBounds.first), log(3.*M.at(i) /(4.*M_PI*VirialOverdensity*DMDensity(z.at(l))))/3. , 1e-4);
				std::cout << SDSBSpline->Val(i, j, l) << std::endl;
			}
		}
	}
	
	SourceDensitySubhaloBoostFT = [SDSBSpline] ( const double k, const double M, const double z)
									{  return SDSBSpline->Eval(M, k, z); };
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
												//std::cout << exp(args[0]) << '\t' << HaloMassFunction(exp(args[0]), params[0]) << '\t' << (1.+ SubhaloBoostFactor(exp(args[0]))) << '\t' << NFW << std::endl;
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
