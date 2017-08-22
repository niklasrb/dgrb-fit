#ifndef BENCHMARKING_CPP
#define BENCHMARKING_CPP

#include "Benchmarking.h"
/// Prepare the Grids for each AstrophysicalSource and then calculate for each individually
void Benchmark::calculateIntensityAndAPSForAstrophysicalSources(std::vector<std::shared_ptr<AstrophysicalSource> > sources)
{
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
		//if(m_log) std::cout << "Starting to calculateIntensityAndAPS for " << source->Name << " at " << t << std::endl;
		calculateIntensityAndAPS(source.get(), EGrid, SGrid, zGrid, GammaGrid);
		//if(m_log) std::cout << std::chrono::duration_cast<std::chrono::millseconds>(std::chrono::system_clock::now()-t).count() << " seconds elapsed" << std::endl; 
	}
	if(m_plot) 
	{
		dIdzCanvas->SetyLimits(1e-8, 1e3);
		(*dIdzCanvas)().BuildLegend();
		(*dIdzCanvas)().SaveAs("dIdz.jpg");
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
	std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline = ObtainSoverLMapping(source, zGrid, GammaGrid);
	if(m_log) SoverLSpline->print();
	if(m_log) std::cout << "Going into ObtainEffectiveEnergySpectrum" << std::endl;
	std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline = ObtainEffectiveEnergySpectrum(source, EGrid, zGrid, GammaGrid, SoverLSpline);
	if(m_log) dNdESpline->print();
	if(m_log) std::cout<< "Going into ObtainFluxThreshold" << std::endl; //std::cin >> dummy;
	std::shared_ptr<gsl2DInterpolationWrapper> S_tIntensitySpline = ObtainFluxThreshold(source, IntensityBins, GammaGrid, SoverLSpline, dNdESpline);
	if(m_log) S_tIntensitySpline->print();
	std::shared_ptr<gsl2DInterpolationWrapper> S_tAPSSpline = ObtainFluxThreshold(source, APSBins, GammaGrid, SoverLSpline, dNdESpline);
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
							return exp(args[0]*(1.+params[1])) * dNdSSpline->Eval( exp(args[0]), params[0]); },
						SGrid.at(0), SGrid.at(SGrid.size()-1) , 2);  			// check range
	
	TF1* GammaIntegrand = new TF1((std::string("Integrated S^alpha * dN/dS for ") + source->Name).c_str(),
										[ SbdNdS, S_tIntensitySpline, &SGrid] (double* args, double* params) // args[0]: Gamma   params[0]: EBin 
										{	SbdNdS->SetParameter(0, args[0]);
											//std::cout << "in Integrated S^alpha * dN/dS: Gamma: " << args[0] << "  S_t = " << S_tSpline->Eval(params[0], args[0]) << std::endl;
											if(S_tIntensitySpline->Eval(params[0], args[0]) <= SGrid.at(0)) return 0.;
											return SbdNdS->Integral(log(SGrid.at(0)), log(S_tIntensitySpline->Eval(params[0], args[0])), 1e-4); },
											source->GammaBounds.first, source->GammaBounds.second, 1);
	//double** C_p = new double*[EBins.size()];
	//for(unsigned int i = 0; i < EBins.size(); i++) C_p[i] = new double[S_t_1.size()];
	
	// first do Intensity
	SbdNdS->SetParameter(1, 1.);  // b =1
	
	for(unsigned int i = 0; i < IntensityBins.size(); i++)
	{
		
		GammaIntegrand->SetParameter(0, (IntensityBins[i].first + IntensityBins[i].second)/2.);
		
		if(GammaGrid.size() == 2) // No Gamma dependency
		{
			const double S = S_tIntensitySpline->Eval((IntensityBins[i].first + IntensityBins[i].second)/2., GammaGrid[0]);
			source->Intensity.push_back(std::pair<Bounds, double>(IntensityBins.at(i), S <= SGrid.at(0)? 0 : SbdNdS->Integral(log(SGrid.at(0)), log(S), 1e-4)  ));
		}
		else
		{
			source->Intensity.push_back(std::pair<Bounds, double>(IntensityBins.at(i), GammaIntegrand->Integral(source->GammaBounds.first, source->GammaBounds.second, 1e-4)));
		}
	}
	
	delete GammaIntegrand;
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

	if(m_plot) plotIntermediates(source, dNdESpline);
	
	/*std::cout << "Ebin\t Intensity\t APS" << std::endl; 
	//for(unsigned int i = 0; i < EBins.size(); i++) std::cout << std::get<1>(EBins[i]) << '\t' << Intensity[i] << '\t' << APS[i] << std::endl;
	
	//std::vector<double> EBinMid; EBinMid.resize(EBins.size()); for(unsigned int i = 0; i < EBins.size(); i++) EBinMid[i] = std::get<1>(EBins[i]);
	/// Now interpolate
	//auto IntensitySpline = std::make_shared<TSpline3>((std::string("Intensity function for ") + source->Name).c_str(),
	//													(double*)EBinMid.data(), (double*)Intensity.data(), EBins.size());
	//auto APSSpline = std::make_shared<TSpline3>((std::string("APS for ") + source->Name).c_str(),
	//													(double*)EBinMid.data(), (double*)APS.data(), EBins.size());													
	//source->Intensity = [IntensitySpline] (const double E) { return IntensitySpline->Eval(E); };
	//source->APS = [APSSpline] (const double E) { return APSSpline->Eval(E); };*/
	
}

/// Calculates dN/dS and saves it in source->dNoverdS(double* args)  args[0]: S, args[1]: Gamma
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtaindNoverdS(AstrophysicalSource* source, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline)
{
	auto dNdSSpline = std::make_shared<gsl2DInterpolationWrapper>(SGrid.data(), SGrid.size(), GammaGrid.data(), GammaGrid.size());
	
	Bounds zBounds_local; zBounds_local.first = std::max(source->zBounds.first, zBounds_global.first);
	zBounds_local.second = std::min(source->zBounds.second, zBounds_global.second);
	
	const double deltaS_min = std::max(SGrid.at(0)*0.1, 1e-70); 	// since deltaS = 0.1*S  we have to troubleshoot for S=0
	
	// Since the integral borders of the inner integral depend on the outer integral, we have to do it like this
	
	TF1* rohdVdz = new TF1((std::string("Integrand rho dV/dz for ") + source->Name).c_str(),   
							[source, this] (double *args, double* params) // args[0]: log(L)  params[0]:z   params[1]: Gamma   
							{ 	//std::cout << source->RescaledLuminosityFunction(args[0], params[0], params[1]) << '\t' << CM->ComovingVolume(params[0]) << std::endl;
								//if(isnan(source->RescaledLuminosityFunction(args[0], params[0], params[1]))) std::cout << "nan for L = " << args[0] << "  z = " << params[0] << "  Gamma = " << params[1] << std::endl;
								return  exp(args[0])*source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1])*CM->ComovingVolume(params[0]); },
							log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 2);
							
	TF1* IntegratedrohdVdzOverL = new TF1((std::string("Integrated(rho dV/dz) over L for ") + source->Name).c_str(),   
										[source, rohdVdz, SoverLSpline, deltaS_min] (double *args, double* params) // args[0]: log(z)  params[0]: S  params[1]: Gamma 
										{ 	double SoverL = SoverLSpline->Eval(exp(args[0]), params[1]);
											rohdVdz->SetParameters(exp(args[0])/*z*/, params[1]/*Gamma*/);
											//std::cout << "S: " << params[0] <<"SOverL: " << SoverL << "   Int bounds: " << params[0] / SoverL << " - " << std::max(1.1*params[0], params[0] + deltaS_min)/ SoverL << std::endl;
											return  exp(args[0])*rohdVdz->Integral(log(params[0] / SoverL) , log(std::max(1.1*params[0], params[0] + deltaS_min)/ SoverL), 1e-4); }, // Integrate from log(L(S)) to log(L(S+deltaS)) 
										zBounds_local.first, zBounds_local.second, 2);
	for(unsigned int i = 0; i < SGrid.size() ; i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			IntegratedrohdVdzOverL->SetParameters(SGrid[i], GammaGrid[j]);
			dNdSSpline->Val(i, j)  = IntegratedrohdVdzOverL->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-4)/(std::max(0.1*SGrid[i], deltaS_min), 1e-4);
			//std::cout << "(" << i << ", " << j << ", " << dNdS[i][j] <<")" << std::endl;
			// delta_S = 0.1*S
		}
	}
	delete IntegratedrohdVdzOverL;
	delete rohdVdz;
	
	dNdSSpline->Initialize();	

	//source->dNoverdS = [dNdSSpline] (const double S, const double Gamma, const double* other) // Save interpolated function in source
	//													{ 	return dNdSSpline->Eval(S, Gamma); };
	return dNdSSpline;
}										

/// Obtain mapping between flux S and Luminosity L for different redshifts and - if applicable - different photon indeces
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainSoverLMapping(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid)
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
								{ return args[0]*source->EnergySpectrum(args[0], params[0], params[1]) / source->kCorrection(args[0], params[0], params[1]) ; },
								/*Emin*/ 0.1_GeV, /*Emax*/ 100._GeV, /*npar*/ 2);
	
	
	// Integrate
	for(unsigned int i = 0; i < zGrid.size(); i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			SIntegrand->SetParameters(zGrid[i], GammaGrid[j]);
			LIntegrand->SetParameters(zGrid[i], GammaGrid[j]);
			SoverLSpline->Val(i, j) = 
					SIntegrand->Integral(0.1*GeV, 100*GeV, 1e-4) / (LIntegrand->Integral(0.1*GeV, 100*GeV, 1e-4)*4*M_PI* pow(CM->ComovingDistance(zGrid[i]),2.));
			// S / (L * 4*pi*d_L^2 )
		}
	}
	delete SIntegrand; delete LIntegrand;  				// Clean up	
	
	// Interpolate S/L over the Grid
	SoverLSpline->Initialize();

	return SoverLSpline;
}

std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<Bounds>& EBins, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline)
{	
	// Calculate flux threshold
	std::vector<double> EBinsMid; EBinsMid.resize(EBins.size());
	for(unsigned int i = 0; i < EBinsMid.size(); i++) EBinsMid[i] = (EBins[i].first + EBins[i].second)/2.;
	
	auto S_tSpline = std::make_shared<gsl2DInterpolationWrapper>(EBinsMid.data(), EBins.size(), GammaGrid.data(), GammaGrid.size());
	
	TF1* dNdE = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(),
							[dNdESpline, source] (double *args, double* params) // args[0]: Energy   params[0]: Gamma
							{ return  std::max(0., dNdESpline->Eval(args[0], params[0])); },
							EBins[0].first, EBins[EBins.size()-1].second, 1);
	
	double denominator = 0;
	for(unsigned int i = 0; i < GammaGrid.size(); i++)
	{
		dNdE->SetParameter(0, GammaGrid[i]);
		denominator = dNdE->Integral(1._GeV, 100._GeV); 
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
									//std::cout << "1D Integrand: " << exp(args[0]) << '\t' << S << '\t' << CM->ComovingVolume(params[0]) << '\t' << source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1]) << '\t' << (1.-DT->DetectionEfficiency(S)) << std::endl;
								  return exp(args[0])*S*CM->ComovingVolume(params[0])*source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1])*(1.-DT->DetectionEfficiency(S)); },
								 log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), /*npar*/ 2);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			Integrand->SetParameters(zGrid[i], GammaGrid[0]);
			dIoverdz.at(i) = Integrand->Integral(log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 1e-4);
		}		
		delete Integrand;												
	}
	else
	{
		TF2* Integrand = new TF2("2D Integrand for dI/dz",     
								[source, SoverLSpline, this] (double* args, double* params) // args[0]: log(L)  args[1]: Gamma   params[0]: z 
								{ double S = exp(args[0])* SoverLSpline->Eval(params[0], args[1]);
									double val =  exp(args[0]) * S * CM->ComovingVolume(params[0]) * source->RescaledLuminosityFunction(exp(args[0]), params[0], args[1]) * (1.-DT->DetectionEfficiency(S));
									//std::cout << exp(args[0]) << '\t' << S << '\t' << CM->ComovingVolume(params[0]) << '\t' << source->RescaledLuminosityFunction(expf(args[0]), params[0], args[1]) << '\t' << (1.-D->DetectionEfficiency(S)) << std::endl;
									 return val; },
								 log(LuminosityBounds_global.first), log(LuminosityBounds_global.second),
								 source->GammaBounds.first, source->GammaBounds.second, /*npar*/ 1);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			//std::cout << zGrid[i] << " - " << CM->ComovingVolume(zGrid[i]) << std::endl;
			Integrand->SetParameters(zGrid[i], 0.);
			dIoverdz.at(i) = Integrand->Integral(log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 
												source->GammaBounds.first, source->GammaBounds.second , 1e-4);
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
		auto g = new TGraph(dIdzSpline->MakeGraph()); g->SetLineColor(std::rand() % 20 + 1);
		dIdzCanvas->AddGraph(g, source->Name, "L");
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
				dNdESpline->Val(j, i) = dIdzdNdE->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-4)/denominator;;
			}
		}
	}
	delete dIdzdNdE;
	delete dIdz;

	dNdESpline->Initialize();
	if(m_log) dNdESpline->print(); // check
	
	return dNdESpline;
}

void Benchmark::plotIntermediates(AstrophysicalSource* source, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline)
{
	return;
}

/// Darkmatter

void Benchmark::calculateIntensityAndAutocorrelationForDM(std::vector<std::shared_ptr<DarkMatter> > DM)
{
	assert(zGridLen >=1);
	std::vector<double> zGrid; zGrid.resize(zGridLen);
	std::vector<double> kGrid; kGrid.resize(kGridLen);
	std::vector<int> Multipoles; Multipoles.resize(5);
	
	for(unsigned int i = 0; i < zGridLen; i++) zGrid.at(i) = exp(log(zBounds_global.first) + double(i)* (log(zBounds_global.second) - log(zBounds_global.first))/(zGridLen-1.));
	for(unsigned int i = 0; i < kGridLen; i++) kGrid.at(i) = exp(log(kBounds_global.first) + double(i)* (log(kBounds_global.second) - log(kBounds_global.first))/(kGridLen-1.));
	for(unsigned int i = 0; i < Multipoles.size(); i++) Multipoles.at(i) = i;
	
	std::cout << zGrid[0] << zGrid[zGridLen-1] << std::endl;
	
	for(unsigned int i = 0; i < DM.size(); i++)
	{
		calculateIntensityForDM(DM.at(i), IntensityBins);
		calculateAPSForDM(DM.at(i), APSBins, zGrid, kGrid, Multipoles);
	}
}

void Benchmark::calculateIntensityForDM(std::shared_ptr<DarkMatter> DM, const std::vector<Bounds>& EBins)
{
	/// Calculate Intensity by integrating window function
	TF2* wf = new TF2((std::string("Window function for") + DM->Name).c_str(),
						[DM] (double* args, double* params) // args[0]: E  args[1]: z
						{	return DM->WindowFunction(args[0], args[1]); 	},
						EBounds_global.first, EBounds_global.second, zBounds_global.first, zBounds_global.second, 0);
	//std::vector<double> EMid; EMid.resize(EBins.size());
	for(unsigned int i = 0; i < EBins.size(); i++)
	{
		//EMid.at(i) = std::get<1>(EBins.at(i));
		DM->Intensity.push_back(std::pair<Bounds, double>(EBins.at(i), wf->Integral(EBins.at(i).first, EBins.at(i).second, zBounds_global.first, zBounds_global.second, 1e-4)));
	}
	delete wf;
}

void Benchmark::calculateAPSForDM(std::shared_ptr<DarkMatter> DM, const std::vector<Bounds>& EBins, const std::vector<double>& zGrid, const std::vector<double>& kGrid, const std::vector<int>& Multipoles)
{
	/// Calculate 3D Power Spectrum and then APS
	TF1* P1HaloIntegrand = new TF1((std::string("The Integrand dn/dm * sourcedensityFT^2 for the 1 Halo term for ") + DM->Name).c_str(),
									[DM, this] (double* args, double* params) // args[0]: log(M)  params[0]: k  params[1]: z
									{ //std::cout << "P1: " << exp(args[0]) << '\t' << params[0] << '\t' << HM->HaloMassFunction(exp(args[0]), params[1]) << '\t' << DM->SourceDensityFT(params[0], exp(args[0]), params[1]) << std::endl;
										return exp(args[0])*HM->HaloMassFunction(exp(args[0]), params[1]) * pow(DM->SourceDensityFT(params[0], exp(args[0]), params[1]), 2); },
									log(HM->MBounds.first), log(HM->MBounds.second), 2);
	
	TF1* P2HaloIntegrand = new TF1((std::string("The Integrand dn/dm * LinearHaloBias *sourcedensityFT for the 2 Halo term for ") + DM->Name).c_str(),
									[DM, this] (double* args, double* params) // args[0]: log(M)  params[0]: k  params[1]: z
									{ //std::cout << "P2: " << exp(args[0]) << '\t' << params[0] << '\t' << HM->HaloMassFunction(exp(args[0]), params[1]) << '\t' << HM->LinearHaloBias(exp(args[0]), params[1]) << '\t' << DM->SourceDensityFT(params[0], exp(args[0]), params[1]) << std::endl;
										return exp(args[0])*HM->HaloMassFunction(exp(args[0]), params[1]) * HM->LinearHaloBias(exp(args[0]), params[1]) * DM->SourceDensityFT(params[0], exp(args[0]), params[1]); },
									log(HM->MBounds.first), log(HM->MBounds.second), 2);
	
	double** _3DPowerSpectrum/*[kGrid.size()][zGrid.size()];*/ = new double*[kGrid.size()]; for(unsigned int i = 0; i < kGrid.size(); i++) _3DPowerSpectrum[i] = new double[zGrid.size()];
	for(unsigned int i = 0; i < kGrid.size(); i++)
	{
		for(unsigned int j = 0; j < zGrid.size(); j++)
		{
			P1HaloIntegrand->SetParameters(kGrid.at(i), zGrid.at(j));
			P2HaloIntegrand->SetParameters(kGrid.at(i), zGrid.at(j));
			_3DPowerSpectrum[i][j] = P1HaloIntegrand->Integral(log(HM->MBounds.first), log(HM->MBounds.second), 1e-4)   
										+ pow(P2HaloIntegrand->Integral(log(HM->MBounds.first), log(HM->MBounds.second), 1e-4), 2)*(*(HM->Plin))(kGrid.at(i), zGrid.at(j));
			
			//std::cout << '(' << i << ", " << j << "): " << _3DPowerSpectrum[i][j] << std::endl;
		}
	} 
	delete P1HaloIntegrand; delete P2HaloIntegrand;
	
	auto _3DPowerSpectrumSpline = std::make_shared<gsl2DInterpolationWrapper>(kGrid.data(), kGrid.size(), zGrid.data(), zGrid.size(),(const double**) _3DPowerSpectrum);
	
	for(unsigned int i = 0; i < kGrid.size(); i++) delete[] _3DPowerSpectrum[i];
	delete []_3DPowerSpectrum; 
	_3DPowerSpectrumSpline->print();
	
	
	TF3* APSCrossIntegrand = new TF3((std::string("Integrand WF*WF * P_ij / chi^2 for the APS of") + DM->Name).c_str(),
								[DM, _3DPowerSpectrumSpline, this] (double* args, double *params) // args[0]:log(z) args[1]:E1  args[2]:E2   params[0]: multipole
								{ 	double chi = CM->ComovingDistance(exp(args[0]));
									return exp(args[0])*c_0/(pow(chi, 2.)*CM->HubbleRate(exp(args[0])))*DM->WindowFunction(args[1], exp(args[0]))*DM->WindowFunction(args[2], exp(args[0])) * _3DPowerSpectrumSpline->Eval(params[1]/chi, exp(args[0])); },
								log(zBounds_global.first), log(zBounds_global.second), EBounds_global.first, EBounds_global.second, EBounds_global.first, EBounds_global.second, 1);
	
	auto APS = std::make_shared<AngularPowerSpectrum<double> >(EBins.size(), EBins.size(), Multipoles.size());
	
	for(unsigned int i = 0; i < EBins.size(); i++)
	{
		for(unsigned int j = 0; j < EBins.size(); j++)
		{
			if(i != j)
			{
				for(unsigned int k = 0; k < Multipoles.size(); k++)
				{
					APSCrossIntegrand->SetParameter(0, Multipoles.at(k));
					(*APS)(i, j, k) = APSCrossIntegrand->Integral(log(zBounds_global.first), log(zBounds_global.second), EBins[i].first, EBins[i].second,
																EBins[j].first, EBins[j].second, 1e-4);
					//std::cout << "( " << i << ", " << j << ", " << k << "): " << (*APS)(i,j,k) << std::endl;
				}
			}
		}
	}
	delete APSCrossIntegrand;
	TF2* APSAutoIntegrand = new TF2((std::string("Integrand WF^2 * P_ij / chi^2 for the APS of") + DM->Name).c_str(),
								[DM, _3DPowerSpectrumSpline, this] (double* args, double *params) // args[0]:log(z) args[1]:E   params[0]: multipole
								{ 	double chi = CM->ComovingDistance(exp(args[0]));
									return exp(args[0])*c_0/(pow(chi, 2.)*CM->HubbleRate(exp(args[0])))*pow(DM->WindowFunction(args[1], exp(args[0])), 2)* _3DPowerSpectrumSpline->Eval(params[1]/chi, exp(args[0])); },
								log(zBounds_global.first), log(zBounds_global.second), EBins[0].first, EBins[EBins.size()-1].second, 1);
	for(unsigned int i = 0; i < EBins.size(); i++)
	{
		for(unsigned int k = 0; k < Multipoles.size(); k++)
		{
			APSAutoIntegrand->SetParameter(0, Multipoles.at(k));
			(*APS)(i, i, k) = APSAutoIntegrand->Integral(log(zBounds_global.first), log(zBounds_global.second), EBins[i].first, EBins[i].second, 1e-4);
			//std::cout << "( " << i << ", " << k << "): " << (*APS)(i,i,k) << std::endl;
		}
	}
	DM->APS = APS;
}

#endif
