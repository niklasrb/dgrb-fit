#ifndef BENCHMARKING_CPP
#define BENCHMARKING_CPP

#include "Benchmarking.h"
/// Prepare the Grids for each AstrophysicalSource and then calculate for each individually
void Benchmark::calculateIntensityAndAPSForAstrophysicalSources(std::vector<std::shared_ptr<AstrophysicalSource> > sources)
{
	if(sources.size() == 0) return;
	assert(zGridLen >= 2);	
	assert(GammaGridLen >= 3);		// To keep GammaGrid.size() == 2 a special case for no gamma dependency	
	assert(SGridLen >= 2);
	
	std::vector<double> SGrid; SGrid.resize(SGridLen);
	for(unsigned int i = 0; i < SGrid.size(); i++) SGrid[i] = exp( log(SBounds_global.first) + i*(log(SBounds_global.second) - log(SBounds_global.first))/(SGrid.size()-1.));
	
	std::vector<double> EGrid; EGrid.resize(EGridLen);
	for(unsigned int i = 0; i < EGrid.size(); i++) EGrid[i] = exp( log(EBounds_global.first) + i*(log(EBounds_global.second) - log(EBounds_global.first))/(EGrid.size()-1.));
	
	std::vector<double> zGrid; zGrid.resize(zGridLen);
	
	std::vector<double> GammaGrid;
	for(unsigned int i = 0; i < sources.size(); i++)
	{
		std::shared_ptr<AstrophysicalSource> source = sources.at(i);
		Bounds zBounds_local; zBounds_local.first = std::max(zBounds_global.first, source->zBounds.first);
		zBounds_local.second = std::min(zBounds_global.second, source->zBounds.second);
		for(unsigned int i = 0; i < zGridLen; i++) 		// logarithmic spacing in redshift
			zGrid[i] = exp(log(zBounds_local.first) + i*(log(zBounds_local.second) - log(zBounds_local.first))/(zGridLen-1));
		if(source->GammaBounds.first >= source->GammaBounds.second) 
		{
			GammaGrid.resize(2);	// we need this for gsl to interpolate in 2D
			GammaGrid[0] = source->GammaBounds.first;
			GammaGrid[1] = GammaGrid[0]+1e-10; 	// arbitrary higher value
		}
		else
		{
			GammaGrid.resize(GammaGridLen);
			for(unsigned int i = 0; i < GammaGridLen; i++) GammaGrid[i] = source->GammaBounds.first + i*(source->GammaBounds.second - source->GammaBounds.first)/(GammaGridLen-1);
		}
		//auto t = std::chrono::system_clock::now();
		if(m_log) std::cout << "Starting to calculateIntensityAndAPS for " << source->Name /*<< " at " << t*/ << std::endl;
		calculateIntensityAndAPS(source.get(), EGrid, SGrid, zGrid, GammaGrid);
		//if(m_log) std::cout << std::chrono::duration_cast<std::chrono::millseconds>(std::chrono::system_clock::now()-t).count() << " seconds elapsed" << std::endl; 
	}
}

/// Calculate the Intensity and Autocorrelation for an individual Astrophysical Source
void Benchmark::calculateIntensityAndAPS(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& SGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid)
{
	assert(source != NULL);
	
	//std::vector<double> SGrid = { 1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4};       // <- Implement somewhere else
	//std::vector<double> S_t_1 = {1e-10, 2e-10, 3e-10, 4e-10, 5e-10, 6e-10, 7e-10, 8e-10, 9e-10, 10e-10}; 
	//double S_t_1GeV = 4.5e-10;
	
	assert(zGrid.size() >= 2);
	assert(GammaGrid.size() >= 2);
	assert(SGrid.size() >= 2);	  
	
	//std::string dummy;
	
	if(m_log) std::cout<< "Going into ObtainSoverLMapping" << std::endl; //std::cin >> dummy;
	std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline = ObtainSoverLMapping(source, Bounds(0.1_GeV, 100._GeV), zGrid, GammaGrid);
	if(m_log) SoverLSpline->print();
	if(m_log) std::cout << "Going into ObtainEffectiveEnergySpectrum" << std::endl;
	std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline = ObtainEffectiveEnergySpectrum(source, EGrid, zGrid, GammaGrid, SoverLSpline);
	if(m_log) dNdESpline->print();
	if(m_log) std::cout<< "Going into ObtainFluxThreshold" << std::endl; //std::cin >> dummy;
	std::shared_ptr<gsl2DInterpolationWrapper> S_tIntensitySpline = ObtainFluxThreshold(source, IntensityBins, GammaGrid, dNdESpline);
	if(m_log) S_tIntensitySpline->print();
	std::shared_ptr<gsl2DInterpolationWrapper> S_tAPSSpline = ObtainFluxThreshold(source, APSBins, GammaGrid, dNdESpline);
	if(m_log) { std::cout << std::endl; S_tAPSSpline->print(); }
	if(m_log) std::cout<< "Going into ObtaindNoverdS" << std::endl; //std::cin >> dummy;
	std::shared_ptr<gsl2DInterpolationWrapper> dNdSSpline = ObtaindNoverdS(source, SGrid, GammaGrid, SoverLSpline);
	if(m_log) dNdSSpline->print();
	if(m_log)std::cout<< "Calculate Intensity and autocorrelation for different Energy bins" << std::endl; //std::cin >> dummy;
	/// Calculate Intensity and autocorrelation for different Energy bins
	
	
	TF1* SbdNdS = new TF1((std::string("Integrand S^alpha * dN/dS for ") + source->Name).c_str(),
						[dNdSSpline] (double* args, double* params) // args[0]: log(S)   params[0]: Gamma  params[1]: alpha
						{ 	//std::cout << "in S^alpha * dN/dS  S = " << args[0] << '\t' << powf(args[0], params[1]) <<  '\t' << source->dNoverdS( args[0], params[0], NULL) <<  std::endl;
							//if(powf(args[0], params[1]
							return exp(args[0]*(1.+params[1])) * std::max(0., dNdSSpline->Eval( exp(args[0]), params[0])); },
						SBounds_global.first, SBounds_global.second, 2);  			// check range
	
	TF1* GammaIntegrand = new TF1((std::string("Integrated S^alpha * dN/dS for ") + source->Name).c_str(),
										[ SbdNdS, S_tIntensitySpline, &SGrid] (double* args, double* params) // args[0]: Gamma   params[0]: EBin 
										{	SbdNdS->SetParameter(0, args[0]);
											//std::cout << "in Integrated S^alpha * dN/dS: Gamma: " << args[0] << "  S_t = " << S_tSpline->Eval(params[0], args[0]) << std::endl;
											if(S_tIntensitySpline->Eval(params[0], args[0]) <= SGrid.at(0)) return 0.;
											return SbdNdS->Integral(log(SGrid.at(0)), log(S_tIntensitySpline->Eval(params[0], args[0])), 1e-4); },
											source->GammaBounds.first, source->GammaBounds.second, 1);
	
	// first do Intensity
	SbdNdS->SetParameter(1, 1.);  // b =1
	
	for(unsigned int i = 0; i < IntensityBins.size(); i++)
	{
		
		GammaIntegrand->SetParameter(0, (IntensityBins[i].first + IntensityBins[i].second)/2.);
		
		if(GammaGrid.size() == 2) // No Gamma dependency
		{
			const double S = S_tIntensitySpline->Eval((IntensityBins[i].first + IntensityBins[i].second)/2., GammaGrid[0]);
			source->Intensity.push_back(S <= SGrid.at(0)? 0 : SbdNdS->Integral(log(SGrid.at(0)), log(S), 1e-4)  );
		}
		else
		{
			source->Intensity.push_back( GammaIntegrand->Integral(source->GammaBounds.first, source->GammaBounds.second, 1e-4));
		}
	}
	
	delete GammaIntegrand;
	
	if(m_log)
	{ std::cout << "Intensity: "; for(unsigned int i = 0; i < IntensityBins.size(); i++) std::cout << "(" << IntensityBins[i].first << ", " << IntensityBins[i].second << "): " << source->Intensity[i] << (i< IntensityBins.size()-1 ? '\t' : '\n'); }
	// now APS
	
	GammaIntegrand = new TF1((std::string("Integrated S^alpha * dN/dS for ") + source->Name).c_str(),
										[ SbdNdS, S_tAPSSpline, &SGrid] (double* args, double* params) // args[0]: Gamma   params[0]: EBin params[1]: S_tSplineMultiplier
										{	SbdNdS->SetParameter(0, args[0]);
											//std::cout << "in Integrated S^alpha * dN/dS: Gamma: " << args[0] << "  S_t = " << S_tSpline->Eval(params[0], args[0]) << std::endl;
											if(params[1]*S_tAPSSpline->Eval(params[0], args[0]) <= SGrid.at(0)) return 0.;
											return SbdNdS->Integral(log(SGrid.at(0)), log(params[1]*S_tAPSSpline->Eval(params[0], args[0])), 1e-4); },
											source->GammaBounds.first, source->GammaBounds.second, 2);
	
	source->APS = std::make_shared<AstrophysicalSourceAPS<std::shared_ptr<gsl1DInterpolationWrapper> > >(APSBins);
	
	SbdNdS->SetParameter(1, 2.);  // b =2
	
	for(unsigned int i = 0; i < APSBins.size(); i++)
	{	
		GammaIntegrand->SetParameter(0, (APSBins[i].first + APSBins[i].second)/2.);
		
		double* C_p = new double[DT->S_t_1.size()];  // Variate S_t to test dependence on flux threshold
		for(unsigned int j = 0; j < DT->S_t_1.size(); j++)
		{
			if(GammaGrid.size() == 2)   
			{
				const double S = DT->S_t_1.at(j)/DT->S_t_1GeV*S_tAPSSpline->Eval((APSBins[i].first + APSBins[i].second)/2., GammaGrid[0]);
				C_p[j] = (S <= SGrid.at(0) ? 0 : SbdNdS->Integral(log(SGrid.at(0)), log(S), 1e-4));
			}
			else
			{
				GammaIntegrand->SetParameter(1, DT->S_t_1.at(j)/DT->S_t_1GeV);		
				C_p[j] = GammaIntegrand->Integral(source->GammaBounds.first, source->GammaBounds.second, 1e-4);
			}
		}
		auto C_pSpline = std::make_shared<gsl1DInterpolationWrapper>(DT->S_t_1.data(), DT->S_t_1.size(), C_p); 
		source->APS->at(i, i) = C_pSpline;
		delete C_p;
		
	}
	delete GammaIntegrand;

	delete SbdNdS;
	
	for(unsigned int i = 0; i < APSBins.size(); i++)
	{
		for(unsigned int j = 0; j < i; j++)
		{
			source->APS->at(i, j) = std::make_shared<gsl1DInterpolationWrapper>(*source->APS->at(i, i).get() *  *source->APS->at(j, j).get());
			source->APS->at(j, i) = source->APS->at(i, j); 
		}
	}
	
	if(m_plot)	// plot N(>S)
	{
		std::vector<double> S; S.resize(20);
		for(unsigned int i =0; i < S.size(); i++) S[i] =  exp(log(1e-13) + i*(log(1e-6) - log(1e-13))/(S.size()-1.));
		
		auto dNdS = new TF1("dNdS", [dNdSSpline] (double* args, double* params) 	// args[0]:log(S) params[0]: Gamma
											{ return exp(args[0])*dNdSSpline->Eval(exp(args[0]), params[0]); },
								log(SBounds_global.first), log(SBounds_global.second), 1);
		dNdS->SetParameters((source->GammaBounds.first + source->GammaBounds.second)/2., 0);
		std::vector<double> N; N.resize(S.size());
		for(unsigned int i =0; i < S.size(); i++) N[i] = dNdS->Integral(log(S[i]), log(SBounds_global.second), 1e-3);
		TGraph* g = new TGraph(N.size(), S.data(), N.data()); g->SetName(source->Name.c_str());
		NofSFile->cd();
		g->Write();
	}
	if(m_plot)	// plot Intensity
	{
		std::vector<double> IntBinMid; IntBinMid.resize(IntensityBins.size());
		std::vector<double> ScaledInt; ScaledInt.resize(IntensityBins.size());
		for(unsigned int i = 0; i < IntBinMid.size(); i++) IntBinMid[i] = (IntensityBins[i].first + IntensityBins[i].second)/2.;
		for(unsigned int j = 0; j < ScaledInt.size(); j++) ScaledInt[j] = source->Intensity[j]*pow(IntBinMid[j], 2.)/(IntensityBins[j].second - IntensityBins[j].first); 
		TGraph* g = new TGraph(IntensityBins.size(), IntBinMid.data(), ScaledInt.data()); g->SetName(source->Name.c_str());
		IntensityFile->cd();
		g->Write();
	}
	
}



/// Calculates dN/dS 
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtaindNoverdS(AstrophysicalSource* source, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline)
{
	auto dNdSSpline = std::make_shared<gsl2DInterpolationWrapper>(SGrid.data(), SGrid.size(), GammaGrid.data(), GammaGrid.size());
	
	Bounds zBounds_local; zBounds_local.first = std::max(source->zBounds.first, zBounds_global.first);
	zBounds_local.second = std::min(source->zBounds.second, zBounds_global.second);
	
	// Evaluate inner integral on a 3dim grid
	//auto rohdVdzIntegrated = std::make_shared<Interpolation3DWrapper>(zGrid.data(), zGrid.size(), GammaGrid.data(), GammaGrid.size(), 
	//Interpolation3DWrapper(const double* _x,unsigned int n_x, const double* _y,unsigned int n_y, const double* _z,unsigned int n_z)
	
	// Since the integral borders of the inner integral depend on the outer integral, we have to do it like this
	
	TF1* rohdVdz = new TF1((std::string("Integrand rho dV/dz for ") + source->Name).c_str(),   
							[source, this] (double *args, double* params) // args[0]: log(L)  params[0]:z   params[1]: Gamma   
							{ 	//std::cout << source->RescaledLuminosityFunction(args[0], params[0], params[1]) << '\t' << CM->ComovingVolume(params[0]) << std::endl;
								//if(isnan(source->RescaledLuminosityFunction(args[0], params[0], params[1]))) std::cout << "nan for L = " << args[0] << "  z = " << params[0] << "  Gamma = " << params[1] << std::endl;
								return  exp(args[0])/1._ergpers*source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1])*CM->ComovingVolumeElement(params[0]); },
							log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 2);
							
	TF1* IntegratedrohdVdzOverL = new TF1((std::string("Integrated(rho dV/dz) over L for ") + source->Name).c_str(),   
										[source, rohdVdz, SoverLSpline] (double *args, double* params) // args[0]: log(z)  params[0]: S  params[1]: Gamma 
										{ 	double SoverL = SoverLSpline->Eval(exp(args[0]), params[1]);
											rohdVdz->SetParameters(exp(args[0])/*z*/, params[1]/*Gamma*/);
											//std::cout << "S: " << params[0] <<"SOverL: " << SoverL << "   Int bounds: " << params[0] / SoverL << " - " << std::max(1.1*params[0], params[0] + deltaS_min)/ SoverL << std::endl;
											return  exp(args[0])*rohdVdz->Integral(log(params[0] / SoverL) , log(1.1*params[0]/SoverL), 1e-4); }, // Integrate from log(L(S)) to log(L(S+deltaS)) 
										zBounds_local.first, zBounds_local.second, 2);
	for(unsigned int i = 0; i < SGrid.size() ; i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			IntegratedrohdVdzOverL->SetParameters(SGrid[i], GammaGrid[j]);
			dNdSSpline->Val(i, j)  = IntegratedrohdVdzOverL->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-4)/9.521e48*pow(1._Mpc,2.)/(0.1*SGrid[i]);	// Mpc^2 to cm^2
			//std::cout << "(" << i << ", " << j << ", " << dNdS[i][j] <<")" << std::endl;
		}
	}
	delete IntegratedrohdVdzOverL;
	delete rohdVdz;
	
	dNdSSpline->Initialize();	

	return dNdSSpline;
}										

/// Obtain mapping between flux S and Luminosity L for different redshifts and - if applicable - different photon indeces
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainSoverLMapping(AstrophysicalSource* source, Bounds EnergyBin, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid)
{
	// InterpolationWrapper that will hold the values				
	auto SoverLSpline = std::make_shared<gsl2DInterpolationWrapper>(zGrid.data(), zGrid.size(), GammaGrid.data(), GammaGrid.size()); 
	
	// Define integrands using ROOT and Lambda functions
	TF1* SIntegrand = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(), 
								[source] (double *args, double *params)   //  args[0]: Energy    params[0]: z   params[1]: Gamma
								{ return source->EnergySpectrum(args[0], params[0], params[1]); },
								/*Emin*/ 0.1_GeV, /*Emax*/ 100._GeV, /*npar*/ 2);
					
	TF1* LIntegrand = new TF1((std::string("Integrand E/k dN/dE for ") + source->Name).c_str(),
								[source] (double *args, double *params)   //  args[0]: Energy    params[0]: z   params[1]: Gamma
								{ return args[0]*source->EnergySpectrumOverK(args[0], params[0], params[1]) ; },
								/*Emin*/ EnergyBin.first, /*Emax*/ EnergyBin.second, /*npar*/ 2);
	
	
	// Integrate
	double L=0, S = 0;
	for(unsigned int i = 0; i < zGrid.size(); i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			SIntegrand->SetParameters(zGrid[i], GammaGrid[j]);
			LIntegrand->SetParameters(zGrid[i], GammaGrid[j]);
			L = LIntegrand->Integral(EnergyBin.first, EnergyBin.second, 1e-3)*1.602e-3  *1._ergpers;		// convert from GeV to erg
			S = SIntegrand->Integral(0.1_GeV, 100._GeV, 1e-3)*1._photonspercm2s;
			SoverLSpline->Val(i, j) = ( L <= 0 ?  0. : S / (L*4.*M_PI* pow(CM->ComovingDistance(zGrid[i]),2.)*9.521e48/pow(1._Mpc,2.) ) );	// Mpc^2 to cm^2 	// lol
		}
	}
	delete SIntegrand; delete LIntegrand;  				// Clean up	
	
	// Interpolate S/L over the Grid
	SoverLSpline->Initialize();

	return SoverLSpline;
}

std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<Bounds>& EBins, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline)
{	
	// Calculate flux threshold
	std::vector<double> EBinsMid; EBinsMid.resize(EBins.size());
	for(unsigned int i =0; i < EBins.size(); i++) EBinsMid[i] = (EBins[i].first + EBins[i].second)/2.;
	
	auto S_tSpline = std::make_shared<gsl2DInterpolationWrapper>(EBinsMid.data(), EBins.size(), GammaGrid.data(), GammaGrid.size());
	
	TF1* dNdE = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(),
							[dNdESpline, source] (double *args, double* params) // args[0]: Energy   params[0]: Gamma
							{ return  std::max(0., dNdESpline->Eval(args[0], params[0])); },
							EBins[0].first, EBins[EBins.size()-1].second, 1);
	
	double denominator = 0;
	for(unsigned int i = 0; i < GammaGrid.size(); i++)
	{
		dNdE->SetParameter(0, GammaGrid[i]);
		denominator = dNdE->Integral(1._GeV, 100._GeV, 1e-4); 
		for(unsigned int j = 0; j < EBins.size(); j++)
		{
			if(denominator <= 0) 
				S_tSpline->Val(j, i) =0;
			else 	
				S_tSpline->Val(j, i) = dNdE->Integral(EBins[j].first, EBins[j].second, 1e-4) * DT->S_t_1GeV / denominator;
		} 			// Integrate in specific Energy Bin
		
	}
	delete dNdE;
	S_tSpline->Initialize();

	return S_tSpline;
}

std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainEffectiveEnergySpectrum(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline)
{
	std::vector<double> dIoverdz; 
	dIoverdz.resize(zGrid.size());
	//std::string dummy;
	if(m_log)std::cout<< "Obtaining dI/dz" << std::endl; //std::cin >> dummy;
	/// Obtain dI/dz
	if(GammaGrid.size() == 2)	// no Gamma Dependency
	{
		TF1* Integrand = new TF1("1D Integrand for dI/dz",   
								[source, SoverLSpline, this] (double* args, double* params) // args[0]: log(L)   params[0]: z   params[1]: Gamma, but should be irrelevant
								{ double S = exp(args[0])* SoverLSpline->Eval(params[0], params[1]);
									//std::cout << "1D Integrand: " << exp(args[0]) << '\t' << S << '\t' << CM->ComovingVolumeElement(params[0]) << '\t' << source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1]) << '\t' << (1.-DT->DetectionEfficiency(S)) << std::endl;
								  return exp(args[0])*S*CM->ComovingVolumeElement(params[0])*source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1])*(1.-DT->DetectionEfficiency(S)); },
								 log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), /*npar*/ 2);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			Integrand->SetParameters(zGrid[i], GammaGrid[0]);
			dIoverdz.at(i) = Integrand->Integral(log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 1e-4)/9.521e48*pow(1._Mpc,2.);
		}		
		delete Integrand;												
	}
	else
	{
		TF2* Integrand = new TF2("2D Integrand for dI/dz",     
								[source, SoverLSpline, this] (double* args, double* params) // args[0]: log(L)  args[1]: Gamma   params[0]: z 
								{ double S = exp(args[0])* SoverLSpline->Eval(params[0], args[1]);
									double val =  exp(args[0]) * S * CM->ComovingVolumeElement(params[0]) * source->RescaledLuminosityFunction(exp(args[0]), params[0], args[1]) * (1.-DT->DetectionEfficiency(S));
									//std::cout << exp(args[0]) << '\t' << S << '\t' << CM->ComovingVolumeElement(params[0]) << '\t' << source->RescaledLuminosityFunction(exp(args[0]), params[0], args[1]) << '\t' << (1.-DT->DetectionEfficiency(S)) << std::endl;
									 return val; },
								 log(LuminosityBounds_global.first), log(LuminosityBounds_global.second),
								 source->GammaBounds.first, source->GammaBounds.second, /*npar*/ 1);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			Integrand->SetParameters(zGrid[i], 0.);
			dIoverdz.at(i) = Integrand->Integral(log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 
												source->GammaBounds.first, source->GammaBounds.second , 1e-4)/9.521e48*pow(1._Mpc,2.);
		}		
		delete Integrand;
	}
	
	if(m_log)
	{
		std::cout << "dI/dz = ";
		for(unsigned int i = 0; i < zGrid.size(); i++) std::cout << "(" << zGrid[i] << "," <<dIoverdz[i] << ")" << '\t';
	}
	
	/// Use dI/dz to calculate effective energy spectrum
	auto dIdzSpline = std::make_shared<gsl1DInterpolationWrapper>(zGrid.data(), zGrid.size(), dIoverdz.data(), gsl_interp_linear, 0); 
	
	if(m_plot)
	{
		auto g = new TGraph(dIdzSpline->MakeGraph());  g->SetName(source->Name.c_str());
		dIdzFile->cd();
		g->Write();
	}	
	
	if(m_log) std::cout<< std::endl << "calculate effective energy spectrum" << std::endl; //std::cin >> dummy;
	Bounds zBounds_local; zBounds_local.first = zGrid[0];
	zBounds_local.second = zGrid[zGrid.size()-1];
	
	auto dNdESpline = std::make_shared<gsl2DInterpolationWrapper>(EGrid.data(), EGrid.size(), GammaGrid.data(), GammaGrid.size());
	
	// Integrate it to get dN/dE without z dependence
	
	TF1* dIdzdNdE = new TF1((std::string("Integrand dI/dz * dN/dE for ") + source->Name).c_str(),
							[dIdzSpline, source] (double* args, double* params) // args[0]: log z  params[0]:E  params[1]:Gamma
							{ //std::cout << "z = " << exp(args[0]) << "  E = " << params[0] << "  Gamma = " << params[1]
								//	<< "   dI/dz = " << dIdzSpline->Eval(exp(args[0])) << "  dN/dE = " << source->EnergySpectrum(params[0], args[0], params[1]) << std::endl;
							 return exp(args[0])* dIdzSpline->Eval(exp(args[0])) * source->EnergySpectrum(params[0], exp(args[0]), params[1]); },
							zBounds_local.first, zBounds_local.second, /*npar*/ 2);
							
							
	TF1* dIdz = new TF1((std::string("Integrand dI/dz for ") + source->Name).c_str() ,
							[dIdzSpline, source] (double* args, double* params) // args[0]: log(z) 
							{ return exp(args[0])*dIdzSpline->Eval(exp(args[0])); },
							zBounds_local.first, zBounds_local.second, /*npar*/ 0);
							
	double denominator = dIdz->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-4);
	for(unsigned int i = 0; i < GammaGrid.size() ; i++)
	{
		for(unsigned int j = 0; j < EGrid.size(); j++)
		{
			if(denominator <= 0) dNdESpline->Val(j, i) = 0;
			else
			{
				dIdzdNdE->SetParameters(EGrid[j], GammaGrid[i]);
				dNdESpline->Val(j, i) = dIdzdNdE->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-3)/denominator;;
			}
		}
	}
	delete dIdzdNdE;
	delete dIdz;

	dNdESpline->Initialize();
	return dNdESpline;
}

/// Darkmatter

void Benchmark::calculateIntensityAndAutocorrelationForDM(std::vector<std::shared_ptr<DarkMatter> >& DM, const std::vector<double>& Multipoles)
{
	assert(zGridLen >=1);
	std::vector<double> zGrid; zGrid.resize(zGridLen);
	std::vector<double> kGrid; kGrid.resize(kGridLen);
	
	for(unsigned int i = 0; i < zGridLen; i++) zGrid.at(i) = exp(log(zBounds_global.first) + double(i)* (log(zBounds_global.second) - log(zBounds_global.first))/(zGridLen-1.));
	for(unsigned int i = 0; i < kGridLen; i++) kGrid.at(i) = exp(log(kBounds_global.first) + double(i)* (log(kBounds_global.second) - log(kBounds_global.first))/(kGridLen-1.));
	
	for(unsigned int i = 0; i < DM.size(); i++)
	{
		std::cout << "calculate Intesity for " << DM[i]->Name << std::endl;
		calculateIntensityForDM(DM[i], IntensityBins);
		std::cout << "calculate APS for " << DM[i]->Name << std::endl;
		calculateAPSForDM(DM[i], APSBins, zGrid, kGrid, Multipoles);
	}
}

void Benchmark::calculateIntensityForDM(std::shared_ptr<DarkMatter> DM, const std::vector<Bounds>& EBins)
{
	/// Calculate Intensity by integrating window function
	TF2* wf = new TF2((std::string("Window function for") + DM->Name).c_str(),
						[DM] (double* args, double* params) // args[0]: E  args[1]: log(z)
						{	return exp(args[1])*DM->WindowFunction(args[0], args[1]); 	},
						EBounds_global.first, EBounds_global.second, log(zBounds_global.first), log(zBounds_global.second), 0);
	
	for(unsigned int i = 0; i < EBins.size(); i++)
	{
		DM->Intensity.push_back( wf->Integral(EBins.at(i).first, EBins.at(i).second, log(zBounds_global.first), log(zBounds_global.second), 1e-4));
	}
	delete wf;
	
	if(m_plot)	// plot Intensity
	{
		std::vector<double> IntBinMid; IntBinMid.resize(EBins.size());
		std::vector<double> ScaledInt; ScaledInt.resize(EBins.size());
		for(unsigned int i = 0; i < IntBinMid.size(); i++) IntBinMid[i] = (EBins[i].first + EBins[i].second)/2.;
		for(unsigned int j = 0; j < ScaledInt.size(); j++) ScaledInt[j] = DM->Intensity[j]*pow(IntBinMid[j], 2.)/(EBins[j].second - EBins[j].first); 
		TGraph* g = new TGraph(IntensityBins.size(), IntBinMid.data(), ScaledInt.data()); g->SetName(DM->Name.c_str());
		IntensityFile->cd();
		g->Write();
	}
}

void Benchmark::calculateAPSForDM(std::shared_ptr<DarkMatter>& DM, const std::vector<Bounds>& EBins, const std::vector<double>& zGrid, const std::vector<double>& kGrid, const std::vector<double>& Multipoles)
{
	/// Calculate 3D Power Spectrum and then APS
	TF1* P1HaloIntegrand = new TF1((std::string("The Integrand dn/dm * sourcedensityFT^2 for the 1 Halo term for ") + DM->Name).c_str(),
									[DM, this] (double* args, double* params) // args[0]: log(M)  params[0]: k  params[1]: z
									{ //std::cout << "P1: " << exp(args[0]) << '\t' << params[0] << '\t' << HM->HaloMassFunction(exp(args[0]), params[1]) << '\t' << DM->SourceDensityFT(params[0], exp(args[0]), params[1]) << std::endl;
										return exp(args[0])*HM->HaloMassFunction(exp(args[0]), params[1]) * pow(DM->SourceDensityFT(params[0], exp(args[0]), params[1]), 2.); },
									log(HM->MBounds.first), log(HM->MBounds.second), 2);
	
	TF1* P2HaloIntegrand = new TF1((std::string("The Integrand dn/dm * LinearHaloBias *sourcedensityFT for the 2 Halo term for ") + DM->Name).c_str(),
									[DM, this] (double* args, double* params) // args[0]: log(M)  params[0]: k  params[1]: z
									{ //std::cout << "P2: " << exp(args[0]) << '\t' << params[0] << '\t' << HM->HaloMassFunction(exp(args[0]), params[1]) << '\t' << HM->LinearHaloBias(exp(args[0]), params[1]) << '\t' << DM->SourceDensityFT(params[0], exp(args[0]), params[1]) << std::endl;
										return exp(args[0])*HM->HaloMassFunction(exp(args[0]), params[1]) * HM->LinearHaloBias(exp(args[0]), params[1]) * DM->SourceDensityFT(params[0], exp(args[0]), params[1]); },
									log(HM->MBounds.first), log(HM->MBounds.second), 2);
	
	auto _3DPowerSpectrumSpline = std::make_shared<gsl2DInterpolationWrapper>(kGrid.data(), kGrid.size(), zGrid.data(), zGrid.size());
	
	double p1, p2;
	for(unsigned int i = 0; i < kGrid.size(); i++)
	{
		for(unsigned int j = 0; j < zGrid.size(); j++)
		{
			P1HaloIntegrand->SetParameters(kGrid.at(i), zGrid.at(j));
			P2HaloIntegrand->SetParameters(kGrid.at(i), zGrid.at(j));
			p1 = P1HaloIntegrand->Integral(log(HM->MBounds.first), log(HM->MBounds.second), 1e-3);
			p2 = P2HaloIntegrand->Integral(log(HM->MBounds.first), log(HM->MBounds.second), 1e-3);
			_3DPowerSpectrumSpline->Val(i, j) =  p1 + pow(p2, 2)*(*(HM->Plin))(kGrid.at(i), zGrid.at(j));
			
			if(m_log) std::cout << '(' << i << ", " << j << "): " << _3DPowerSpectrumSpline->Val(i, j) << std::endl;
		}
	} 
	delete P1HaloIntegrand; delete P2HaloIntegrand;
	_3DPowerSpectrumSpline->Initialize();
	if(m_log) {std::cout << "3D Power Spectrum: "; _3DPowerSpectrumSpline->print(); }
	
	if(m_plot)
	{
		std::vector<double> PS; PS.resize(kGrid.size());
		for(unsigned int i = 0; i < kGrid.size(); i++) PS[i] = _3DPowerSpectrumSpline->Eval(kGrid[i]*h, 1);
		auto g = new TGraph(kGrid.size(), kGrid.data(), PS.data()); g->SetName(DM->Name.c_str());
		_3DPSFile->cd();
		g->Write();
	}
	
	DM->APS = std::make_shared<AngularPowerSpectrum<double> >(EBins.size(), EBins.size(), Multipoles.size());
	
	if(m_log) std::cout << "calculate APS crosscorelation" << std::endl;
	TF3* APSCrossIntegrand = new TF3((std::string("Integrand WF*WF * P_ij / chi^2 for the APS of") + DM->Name).c_str(),
								[DM, _3DPowerSpectrumSpline, this] (double* args, double *params) // args[0]:log(z) args[1]:E1  args[2]:E2   params[0]: multipole
								{ 	double chi = CM->ComovingDistance(exp(args[0]));
									return exp(args[0])*c_0/(pow(chi, 2.)*CM->HubbleRate(exp(args[0])))*DM->WindowFunction(args[1], exp(args[0]))*DM->WindowFunction(args[2], exp(args[0])) * _3DPowerSpectrumSpline->Eval(params[0]/chi, exp(args[0])); },
								log(zBounds_global.first), log(zBounds_global.second), EBounds_global.first, EBounds_global.second, EBounds_global.first, EBounds_global.second, 1);
	
	
	for(unsigned int i = 0; i < EBins.size(); i++)
	{
		for(unsigned int j = 0; j < i; j++)
		{
			for(unsigned int k = 0; k < Multipoles.size(); k++)
			{
				APSCrossIntegrand->SetParameter(0, Multipoles.at(k));
				DM->APS->at(i, j, k) = APSCrossIntegrand->Integral(log(zBounds_global.first), log(zBounds_global.second), EBins[i].first, EBins[i].second,
															EBins[j].first, EBins[j].second, 1e-3);
				DM->APS->at(j, i, k) = DM->APS->at(i, j, k);
				if(m_log) std::cout << "( " << i << ", " << j << ", " << k << "): " << DM->APS->at(i,j,k) << std::endl;
			}
		}
	}
	delete APSCrossIntegrand;
	
	if(m_log) std::cout << "calculate APS autocorrelation" << std::endl;
	TF2* APSAutoIntegrand = new TF2((std::string("Integrand WF^2 * P_ij / chi^2 for the APS of") + DM->Name).c_str(),
								[DM, _3DPowerSpectrumSpline, this] (double* args, double *params) // args[0]:log(z) args[1]:E   params[0]: multipole
								{ 	const double chi = CM->ComovingDistance(exp(args[0]));
									//std::cout << exp(args[0]) << "\t" << args[1] << "\t" << chi << '\t' << params[0]/chi << "\t" << DM->WindowFunction(args[1], exp(args[0])) << "\t" << _3DPowerSpectrumSpline->Eval(params[0]/chi, exp(args[0])) << std::endl;
									return exp(args[0])*c_0/(pow(chi, 2.)*CM->HubbleRate(exp(args[0])))*pow(DM->WindowFunction(args[1], exp(args[0])), 2)* _3DPowerSpectrumSpline->Eval(params[0]/chi, exp(args[0])); },
								log(zBounds_global.first), log(zBounds_global.second), EBounds_global.first, EBounds_global.second, 1);
	for(unsigned int i = 0; i < EBins.size(); i++)
	{
		for(unsigned int k = 0; k < Multipoles.size(); k++)
		{
			APSAutoIntegrand->SetParameter(0, Multipoles[k]);
			DM->APS->at(i, i, k) = APSAutoIntegrand->Integral(log(zBounds_global.first), log(zBounds_global.second), EBins[i].first, EBins[i].second, 1e-3);
			if(m_log) std::cout << "( " << i << ", " << k << "): " << DM->APS->at(i,i,k) << std::endl;
		}
	}
	delete APSAutoIntegrand;
}

#endif
