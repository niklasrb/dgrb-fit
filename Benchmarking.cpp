#ifndef BENCHMARKING_CPP
#define BENCHMARKING_CPP

#include "Benchmarking.h"


/// Obtain mapping between flux S and Luminosity L for different redshifts and - if applicable - different photon indeces
void Benchmark::ObtainSoverLMapping(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl_spline2d>& SoverLSpline, std::shared_ptr<gsl_interp_accel>& zAcc, std::shared_ptr<gsl_interp_accel>& GammaAcc)
{
	int GammaGridSize = GammaGrid.size();
	int zGridSize = zGrid.size();
	// Matrix for S/L
	double* SoverL = new double[GammaGridSize*zGridSize]; 							
	
	// Define integrands using ROOT and Lambda functions
	//  args[0] is Energy    params[0]: z   params[1]: Gamma
	TF1* SIntegrand = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(), 
					[source] (double *args, double *params)  { return source->EnergySpectrum(args[0], params[0], params[1]); },
					/*xmin*/ 0.1_GeV, /*xmax*/ 100._GeV, /*npar*/ 2);
	TF1* LIntegrand = new TF1((std::string("Integrand E/k dN/dE for ") + source->Name).c_str(),
					[source] (double *args, double *params) 
					{ return args[0]*source->EnergySpectrum(args[0], params[0], params[1]) / source->kCorrection(args[0], params[0], params[1]) ; },
					/*xmin*/ 0.1_GeV, /*xmax*/ 100._GeV, /*npar*/ 2);
	
	// Set Pointer to spline for calling function to use
	SoverLSpline.reset(gsl_spline2d_alloc(gsl_interp2d_bilinear, zGrid.size(), GammaGrid.size()), [SoverL](gsl_spline2d* s){gsl_spline2d_free(s); delete []SoverL; std::cout << "SoverLSpline deconstructor" << std::endl;}); // Destructor will be called when the shared_ptr reaches the end of its life
	zAcc.reset(gsl_interp_accel_alloc(), [](gsl_interp_accel* acc) { gsl_interp_accel_free(acc); });
	GammaAcc.reset(gsl_interp_accel_alloc(), [](gsl_interp_accel* acc) { gsl_interp_accel_free(acc); });
	
	// Integrate
	for(int i = 0; i < zGridSize; i++)
	{
		//SoverL[i] = new double[GammaGridSize];
		for(int j = 0; j < GammaGridSize; j++)
		{
			SIntegrand->SetParameters(zGrid[i], GammaGrid[j]);
			LIntegrand->SetParameters(zGrid[i], GammaGrid[j]);
			gsl_spline2d_set (SoverLSpline.get(), SoverL, i, j, 
							SIntegrand->Integral(0.1*GeV, 100*GeV, 1e-6) / (LIntegrand->Integral(0.1*GeV, 100*GeV, 1e-6)*4*M_PI* powf(CM->ComovingDistance(zGrid[i]),2)));
			// S / (L * 4*pi*d_L^2 )
		}
	}
	delete SIntegrand; delete LIntegrand;  				// Clean up
	
	std::cout << "SoverL: "; for(int i = 0; i < zGridSize*GammaGridSize; i++) std::cout << SoverL[i] << '\t';
	
	// Interpolate S/L over the Grid
	gsl_spline2d_init(SoverLSpline.get(), zGrid.data(), GammaGrid.data(), SoverL, zGrid.size(), GammaGrid.size());
}

std::vector<double> Benchmark::ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl_spline2d>& SoverLSpline, std::shared_ptr<gsl_interp_accel>& zAcc, std::shared_ptr<gsl_interp_accel>& GammaAcc)
{
	int zGridSize = zGrid.size(); int GammaGridSize = GammaGrid.size();
	int EBinsSize = EBins.size();
	
	std::vector<double> dIoverdz; 
	dIoverdz.resize(zGridSize);
	std::string dummy;
	std::cout<< "Obtaining dI/dz" << std::endl; //std::cin >> dummy;
	/// Obtain dI/dz
	if(GammaGridSize == 1)
	{
		TF1* Integrand = new TF1("Integrand for dI/dz",   
								[source, SoverLSpline, zAcc, GammaAcc, this] (double* args, double* params) // args[0]: L   params[0]: z   params[1]: Gamma, but should be irrelevant
								{ double S = args[0]*gsl_spline2d_eval(SoverLSpline.get(), params[0], params[1], zAcc.get(), GammaAcc.get());
								  return S*CM->ComovingVolume(params[0])*source->RescaledLuminosityFunction(args[0], params[0], params[1])*(1-D->DetectionEfficiency(S)); },
								 LuminosityBounds.first, LuminosityBounds.second, /*npar*/ 2);
		for(int i = 0; i < zGridSize; i++)
		{
			Integrand->SetParameters(zGrid[i], GammaGrid[0]);
			dIoverdz.at(i) = Integrand->Integral(LuminosityBounds.first, LuminosityBounds.second);
		}		
		delete Integrand;												
	}
	else
	{
		TF2* Integrand = new TF2("Integrand for dI/dz",   // args[0]: L  args[1]: Gamma     params[0]: z   
								[source, SoverLSpline, zAcc, GammaAcc, this] (double* args, double* params) 
								{ double S = args[0]*gsl_spline2d_eval(SoverLSpline.get(), params[0], args[1], zAcc.get(), GammaAcc.get());
									double val =  S * CM->ComovingVolume(params[0]) * source->RescaledLuminosityFunction(args[0], params[0], args[1]) * (1.-D->DetectionEfficiency(S));
									//std::cout << S << '\t' << CM->ComovingVolume(params[0]) << '\t' << source->RescaledLuminosityFunction(args[0], params[0], args[1]) << '\t' << (1.-D->DetectionEfficiency(S)) << std::endl;
									 return val; },
								 LuminosityBounds.first, LuminosityBounds.second,
								 source->m_GammaBounds.first, source->m_GammaBounds.second, /*npar*/ 1);
		for(int i = 0; i < zGridSize; i++)
		{
			//std::cout << zGrid[i] << " - " << CM->ComovingVolume(zGrid[i]) << std::endl;
			Integrand->SetParameters(zGrid[i], 0);
			dIoverdz.at(i) = Integrand->Integral(LuminosityBounds.first, LuminosityBounds.second, 
												source->m_GammaBounds.first, source->m_GammaBounds.second, /*epsrel*/ 1e-4)/1e6;
		}		
		delete Integrand;
	}
	
	
	for(int i = 0; i < zGridSize; i++) std::cout << "(" << zGrid[i] << "," <<dIoverdz[i] << ")" << '\t';
	
	std::cout<< "calculate effective energy spectrum" << std::endl; std::cin >> dummy;
	/// Use dI/dz to calculate effective energy spectrum
	TSpline3* dIdzSpline = new TSpline3((std::string("dI/dz Spline for ") + source->Name).c_str(), 
													(double*)zGrid.data(), (double*)dIoverdz.data(), zGridSize); 
	
	if(m_plot)
	{
		dIdzCanvas->cd();
		dIdzSpline->Draw("L");		
	}
	
	std::vector<double> dNoverdE; 
	dNoverdE.resize(EBinsSize*GammaGridSize);
	std::cout << dNoverdE.size();
	// Integrate it to get dN/dE without z dependence
	std::cout << std::endl << "get dN/dE" << std::endl;
	
	TF1* dIdzdNdE = new TF1((std::string("Integrand dI/dz * dN/dE for ") + source->Name).c_str(),
							[dIdzSpline, source] (double* args, double* params) // args[0]: z  params[0]:E  params[1]:Gamma
							{ return dIdzSpline->Eval(args[0]) * source->EnergySpectrum(params[0], args[0], params[1]); },
							source->m_zBounds.first, source->m_zBounds.second, /*npar*/ 2);
	TF1* dIdz = new TF1((std::string("Integrand dI/dz for ") + source->Name).c_str() ,
							[dIdzSpline, source] (double* args, double* params) // args[0]: z 
							{ return dIdzSpline->Eval(args[0]); },
							source->m_zBounds.first, source->m_zBounds.second, /*npar*/ 0);
	double denominator = dIdz->Integral(source->m_zBounds.first, source->m_zBounds.second);
	for(int i = 0; i < GammaGridSize ; i++)
	{
		for(int j = 0; j < EBinsSize; j++)
		{
			dIdzdNdE->SetParameters(std::get<1>(EBins[j]), GammaGrid[i]);
			dNoverdE.at(i+j*GammaGridSize) = dIdzdNdE->Integral(source->m_zBounds.first, source->m_zBounds.second)/denominator;
		}
	}
	delete dIdzdNdE;
	delete dIdz;
	delete dIdzSpline;
	std::cout << "dNoverdE = ";	for(int i = 0; i < EBinsSize*GammaGridSize; i++) std::cout << dNoverdE[i] << '\t'; std::cout << std::endl;
	
	// Calculate flux threshold
	std::vector<double> S_t;
	S_t.resize(EBinsSize*GammaGridSize);
	std::vector<double> EBinsMid; for(int i = 0; i < EBinsSize; i++) EBinsMid.push_back(std::get<1>(EBins[i]));
	
	for(int i = 0; i < GammaGridSize; i++)
	{
		// Since gsl_spline2d does not allow extrapolation for 2D Grids, we have to use 1D inter- and extrapolation for each Gamma
		TSpline3* dNdESpline = new TSpline3((std::string("dN/dE Spline for ") + source->Name).c_str(),
											(double*) EBinsMid.data(), (double*)&dNoverdE[i*GammaGridSize], EBinsSize);
		TF1* dNdE = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(),
							[dNdESpline, source] (double *args, double* params) // args[0]: Energy   
							{ return  dNdESpline->Eval(args[0]); },
							std::get<0>(EBins[0]), std::get<2>(EBins[EBinsSize-1]), 0);
		denominator = dNdE->Integral(1._GeV, /*std::numeric_limits<double>::infinity()*/ 100._GeV); 
		for(int j = 0; j < EBinsSize; j++)
		{
			S_t.at(i+j*GammaGridSize) = dNdE->Integral(std::get<0>(EBins[j]), std::get<2>(EBins[j])) * D->S_t_1GeV() / denominator;
		}
		delete dNdESpline;
		delete dNdE;
	}
	std::cout <<" S_t = ["; for(int i = 0; i < EBinsSize*GammaGridSize; i++) std::cout << S_t[i] << '\t'; std::cout << "]" << std::endl;
	return S_t;
}

void Benchmark::calculateIntensityAndAutocorrelation(int zGridLen, int GammaGridLen)
{
	assert(zGridLen >= 1);
	assert(GammaGridLen >= 1);
	std::vector<double> zGrid; zGrid.resize(zGridLen);
	std::vector<double> GammaGrid;
	for(auto it = AstrophysicalSources.begin(); it != AstrophysicalSources.end(); ++it)
	{
		std::shared_ptr<AstrophysicalSource> source = *it;
		for(int i = 0; i < zGridLen; i++) zGrid[i] = source->m_zBounds.first + i*(source->m_zBounds.second - source->m_zBounds.first)/(zGridLen-1);
		if(source->m_GammaBounds.first == source->m_GammaBounds.second) 
		{
			GammaGrid.resize(1);
			GammaGrid[0] = source->m_GammaBounds.first;
		}
		else
		{
			GammaGrid.resize(GammaGridLen);
			for(int i = 0; i < GammaGridLen; i++) GammaGrid[i] = source->m_GammaBounds.first + i*(source->m_GammaBounds.second - source->m_GammaBounds.first)/(GammaGridLen-1);
		}
		std::time_t t = std::time(nullptr);
		if(m_log) std::cout << "Starting to calculatedNoverdS for " << source->Name << " at " << std::asctime(std::localtime(&t)) << std::endl;
		calculateIntensityAndAutocorrelation(source.get(), zGrid, GammaGrid);
		if(m_log) std::cout << std::time(nullptr) - t << " seconds elapsed" << std::endl; 
	}
	
	if(m_plot)
	{
		dIdzCanvas->SaveAs("dIdz.root");
	}
}


void Benchmark::calculateIntensityAndAutocorrelation(AstrophysicalSource* source,const  std::vector<double>& zGrid, const std::vector<double>& GammaGrid)
{
	assert(source != NULL);
	
	std::vector<double> SGrid;       // <- Implement somewhere else
	int SGridSize = SGrid.size(); 
	
	int zGridSize = zGrid.size(); 
	int GammaGridSize = GammaGrid.size();
	int EBinsSize = EBins.size();		  	// To avoid unnecessary function calls
	
	std::shared_ptr<gsl_spline2d> SoverLSpline;
	std::shared_ptr<gsl_interp_accel> zAcc;
	std::shared_ptr<gsl_interp_accel> GammaAcc;
	
	std::string dummy;
	
	std::cout<< "Going into ObtainSoverLMapping" << std::endl; //std::cin >> dummy;
	ObtainSoverLMapping(source, zGrid, GammaGrid, SoverLSpline, zAcc, GammaAcc);
	std::cout<< "Going into ObtainFluxThreshold" << std::endl; //std::cin >> dummy;
	std::vector<double> S_t = ObtainFluxThreshold(source, zGrid, GammaGrid, SoverLSpline, zAcc, GammaAcc);
	
	std::cout<< "Calculating dN/dS" << std::endl; std::cin >> dummy;
	/// Calculare dN/dS
	
	std::vector<double> dNdS; 
	dNdS.resize(SGridSize*GammaGridSize);
	
	//Since the integral borders of the inner integral depend on the outer integral, we have to do it like this
	
	TF1* rohdVdz = new TF1((std::string("Integrand rho dV/dz for ") + source->Name).c_str(),   
							[source, this] (double *args, double* params) // args[0]: Luminosity  params[0]:z   params[1]: Gamma   
							{ return  source->RescaledLuminosityFunction(args[0], params[0], params[1])*CM->ComovingVolume(params[0]); },
							LuminosityBounds.first, LuminosityBounds.second, 2);
							
	TF1* IntegratedrohdVdzOverL = new TF1((std::string("Integrate(rho dV/dz) over L for ") + source->Name).c_str(),   
										[source, rohdVdz, SoverLSpline, zAcc, GammaAcc] (double *args, double* params) // args[0]: z  params[0]: S  params[1]: Gamma 
										{ 	double SoverL = gsl_spline2d_eval(SoverLSpline.get(), args[0], params[1], zAcc.get(), GammaAcc.get());
											rohdVdz->SetParameters(args[0]/*z*/, params[1]/*Gamma*/);
											return  rohdVdz->Integral(params[0] / SoverL , 1.1*params[0]/ SoverL); }, // Integrate from L(S) to L(S+deltaS) 
										source->m_zBounds.first, source->m_zBounds.second, 2);
	for(int i = 0; i < GammaGridSize; i++)
	{
		for(int j = 0; j < SGridSize; j++)
		{
			IntegratedrohdVdzOverL->SetParameters(SGrid[j], GammaGrid[i]);
			dNdS.at(i + j*GammaGridSize) = IntegratedrohdVdzOverL->Integral(source->m_zBounds.first, source->m_zBounds.second)/(0.1*SGrid[j]);
		}
	}
	delete IntegratedrohdVdzOverL;
	delete rohdVdz;
	
}


#endif
