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
			GammaGrid[1] = source->GammaBounds.first+1; 	// arbitrary higher value
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
	
	assert(zGrid.size() >= 2);
	assert(GammaGrid.size() >= 2);
	assert(SGrid.size() >= 2);	  
	
	//std::string dummy;

	if(m_log) std::cout << "Going into ObtainEffectiveEnergySpectrum" << std::endl;
	std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline = ObtainEffectiveEnergySpectrum(source, EGrid, zGrid, GammaGrid);
	if(m_log) dNdESpline->print();
	if(m_log) std::cout<< "Going into ObtainFluxThreshold" << std::endl; //std::cin >> dummy;
	std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > S_tIntensitySplines = ObtainFluxThreshold(source, IntensityBins, GammaGrid, dNdESpline);
	//if(m_log) S_tIntensitySpline->print();
	std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > S_tAPSSplines = ObtainFluxThreshold(source, APSBins, GammaGrid, dNdESpline);
	//if(m_log) { std::cout << std::endl; S_tAPSSpline->print(); }
	
	/// Calculate Intensity
	if(m_log) std::cout << "calculate intensity for " << source->Name << std::endl;
	for(unsigned int i = 0; i < IntensityBins.size(); i++)
	{
		auto SoverL = ObtainSoverLMapping(source, IntensityBins.at(i), zGrid, GammaGrid);
		
		/// first integration in flux
		auto SrhodVdzLoverSIntegrand = new TF1(("Integrand S*rho*dV/dz *L/S for " + source->Name).c_str(),
										[SoverL, source, this] (double* args, double* params) 	// args[0]: log(S), params[0]: z : params[1]: Gamma
										{	const double L = exp(args[0])/SoverL->Eval(params[0], params[1]); if(L < source->LumBounds.first || L > source->LumBounds.second) return 0.;
											return exp(2*args[0]) * source->RescaledLuminosityFunction( L, params[0], params[1]) 					// undo unit conversion
																	* CM->ComovingVolumeElement(params[0]) /(SoverL->Eval(params[0], params[1])/1.602e-3 *9.521e48/pow(1._Mpc,2.)); },
											log(SBounds_global.first), log(SBounds_global.second), 2);
		
		auto SrhodVdzLoverSIntegratedSpline = std::make_shared<gsl2DInterpolationWrapper>(zGrid.data(), zGrid.size(), GammaGrid.data(), GammaGrid.size());
		
		double S_t = 0;
		for(unsigned int j = 0; j < zGrid.size(); j++)
		{
			for(unsigned int k = 0; k < GammaGrid.size(); k++)
			{
				SrhodVdzLoverSIntegrand->SetParameters(zGrid.at(j), GammaGrid.at(k));
				S_t = S_tIntensitySplines.at(i)->Eval(GammaGrid.at(k));
				if(S_t <= SBounds_global.first) 
					SrhodVdzLoverSIntegratedSpline->Val(j, k) = 0;
				else 
					SrhodVdzLoverSIntegratedSpline->Val(j, k) = SrhodVdzLoverSIntegrand->Integral(log(SBounds_global.first), log(S_t), 1e-3);
			}
		}
		delete SrhodVdzLoverSIntegrand;
		SrhodVdzLoverSIntegratedSpline->Initialize();
		/// Now integrate in redshift
		auto dIdzIntegrand = new TF1(("Integrand \\int dS S*rho*dV/dz*L/S for " + source->Name).c_str(),
							[SrhodVdzLoverSIntegratedSpline] (double* args, double* params) 	// args[0]: log(z)  params[0]: Gamma
							{ return exp(args[0]) *SrhodVdzLoverSIntegratedSpline->Eval(exp(args[0]), params[0]); },
							log(zBounds_global.first), log(zBounds_global.second), 1);
		
		std::vector<double> dIdzIntegrated; dIdzIntegrated.resize(GammaGrid.size());
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			dIdzIntegrand->SetParameters(GammaGrid.at(j), 0.);
			dIdzIntegrated.at(j) = dIdzIntegrand->Integral(log(zGrid.at(0)), log(zGrid.at(zGrid.size()-1)), 1e-3);	// integrate in local zBounds
		}
		delete dIdzIntegrand;
		
		if(GammaGrid.size() == 2)	// no gamma dependency, we're done
		{
			source->Intensity.push_back(dIdzIntegrated.at(0));
		}
		else 	// integrate over Gamma
		{
			auto dIdzIntegratedSpline = std::make_shared<gsl1DInterpolationWrapper>(GammaGrid.data(), GammaGrid.size(), dIdzIntegrated.data());
			auto dIdG = new TF1(("Integrand \\int dz \\int dS S*rho*dV/dz*L/S for " + source->Name).c_str(),
								[dIdzIntegratedSpline] (double* args, double* params)	// args[0]: Gamma
								{ return dIdzIntegratedSpline->Eval(args[0]); },
								source->GammaBounds.first, source->GammaBounds.second, 0);
			source->Intensity.push_back( dIdG->Integral(source->GammaBounds.first, source->GammaBounds.second, 1e-3));
			delete dIdG;
		}
		if(m_log) std::cout << "(" << IntensityBins.at(i).first << ", " << IntensityBins.at(i).second << "): " << source->Intensity.at(i) << std::endl;
	}
	
	/// Calculate APS
	if(m_log) std::cout << "calculate anisotropy for " << source->Name << std::endl;
	source->APS = std::make_shared<AstrophysicalSourceAPS<std::shared_ptr<gsl1DInterpolationWrapper> > >(APSBins);
	for(unsigned int i = 0; i < APSBins.size(); i++)
	{
		auto SoverL = ObtainSoverLMapping(source, APSBins.at(i), zGrid, GammaGrid);
		
		std::vector<double> S_t_1 = DT->S_t_1;	// test dependency on galazy catalog flux threshold 
		
		/// first integration in flux
		auto SrhodVdzLoverSIntegrand = new TF1(("Integrand S^2*rho*dV/dz *L/S for " + source->Name).c_str(),
										[SoverL, source, this] (double* args, double* params) 	// args[0]: log(S), params[0]: z : params[1]: Gamma
										{	const double L = exp(args[0])/SoverL->Eval(params[0], params[1]); if(L < source->LumBounds.first || L > source->LumBounds.second) return 0.;
											return exp(3.*args[0]) * source->RescaledLuminosityFunction( L, params[0], params[1]) 					// undo unit conversion
																	* CM->ComovingVolumeElement(params[0]) /(SoverL->Eval(params[0], params[1])/1.602e-3 *9.521e48/pow(1._Mpc,2.)); },
											log(SBounds_global.first), log(SBounds_global.second), 2);
		
		auto SrhodVdzLoverSIntegratedSpline = std::make_shared<Interpolation3DWrapper>(zGrid.data(), zGrid.size(), GammaGrid.data(), GammaGrid.size(), S_t_1.data(), S_t_1.size());
		
		double S_t = 0;
		for(unsigned int j = 0; j < zGrid.size(); j++)
		{
			for(unsigned int k = 0; k < GammaGrid.size(); k++)
			{
				SrhodVdzLoverSIntegrand->SetParameters(zGrid.at(j), GammaGrid.at(k));
				S_t = S_tIntensitySplines.at(i)->Eval(GammaGrid.at(k));
				 
				for(unsigned int l = 0; l < S_t_1.size(); l++) 
				{
					if(S_t <= SBounds_global.first) 
						SrhodVdzLoverSIntegratedSpline->Val(j, k, l) = 0;
					else
					{
						SrhodVdzLoverSIntegratedSpline->Val(j, k, l) = SrhodVdzLoverSIntegrand->Integral(log(SBounds_global.first), log(S_t/DT->S_t_1GeV * S_t_1.at(l)), 1e-3);
					}														// renormalize S_t_1 and integrate for every grid point
				}
			}
		}
		delete SrhodVdzLoverSIntegrand;

		/// Now integrate in redshift
		auto dCdzIntegrand = new TF1(("Integrand \\int dS S^2*rho*dV/dz*L/S for " + source->Name).c_str(),
							[SrhodVdzLoverSIntegratedSpline] (double* args, double* params) 	// args[0]: log(z)  params[0]: Gamma  params[1]: S_t_1
							{ return exp(args[0]) *SrhodVdzLoverSIntegratedSpline->Eval(exp(args[0]), params[0], params[1]); },
							log(zBounds_global.first), log(zBounds_global.second), 2);
		
		auto dCdzIntegratedSpline = std::make_shared<gsl2DInterpolationWrapper>(GammaGrid.data(), GammaGrid.size(), S_t_1.data(), S_t_1.size());
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			for(unsigned int k = 0; k < S_t_1.size(); k++)
			{
				dCdzIntegrand->SetParameters(GammaGrid.at(j), S_t_1.at(k));
				dCdzIntegratedSpline->Val(j, k) = dCdzIntegrand->Integral(log(zGrid.at(0)), log(zGrid.at(zGrid.size()-1)), 1e-3);	// integrate in local zBounds
			}
		}
		delete dCdzIntegrand;
		dCdzIntegratedSpline->Initialize();
		
		if(GammaGrid.size() == 2)	// no gamma dependency, we're done
		{
			source->APS->at(i, i) = std::make_shared<gsl1DInterpolationWrapper>(dCdzIntegratedSpline->AlongY(GammaGrid.at(0)));
		}
		else 	// integrate over Gamma
		{
			auto dCdG = new TF1(("Integrand \\int dz \\int dS S^2*rho*dV/dz*L/S for " + source->Name).c_str(),
								[dCdzIntegratedSpline] (double* args, double* params)	// args[0]: Gamma  params[0]: S_t_1
								{ return dCdzIntegratedSpline->Eval(args[0], params[0]); },
								source->GammaBounds.first, source->GammaBounds.second, 1);
			std::vector<double> C_p; C_p.resize(S_t_1.size());
			for(unsigned int j = 0; j < S_t_1.size(); j++)
			{
				dCdG->SetParameters(S_t_1.at(j), 0.);
				C_p.at(j) = dCdG->Integral(source->GammaBounds.first, source->GammaBounds.second, 1e-3);
				
			}
			delete dCdG;
			source->APS->at(i, i) = std::make_shared<gsl1DInterpolationWrapper>(S_t_1.data(), S_t_1.size(), C_p.data());
		}
		if(m_log) std::cout << "(" << APSBins.at(i).first << ", " << APSBins.at(i).second << "): " << source->APS->at(i, i)->Eval(DT->S_t_1GeV) << std::endl;
	}
	/// calculate cross corr
	for(unsigned int i = 0; i < APSBins.size(); i++)
	{
		for(unsigned int j = 0; j < i; j++)
		{
			source->APS->at(i, j) = std::make_shared<gsl1DInterpolationWrapper>( sqrt( *(source->APS->at(i,i).get())  *  *(source->APS->at(j,j).get())  ));
			source->APS->at(j, i) = source->APS->at(i, j);
		}
	}
	
	
	/*if(m_log) std::cout<< "Going into ObtaindNoverdS" << std::endl; //std::cin >> dummy;
	std::shared_ptr<gsl2DInterpolationWrapper> dNdSSpline = ObtaindNoverdS(source, Bounds(0.1_GeV, 100._GeV), zGrid, SGrid, GammaGrid);
	//if(m_log) dNdSSpline->print();
	if(m_log) std::cout << "Obtaining dNdS for each energy bin" << std::endl;
	std::vector<std::shared_ptr<gsl2DInterpolationWrapper> > dNdSSplines; dNdSSplines.resize(IntensityBins.size());
	for(unsigned int i = 0; i < IntensityBins.size(); i++) dNdSSplines[i] = ObtaindNoverdS(source, IntensityBins[i], zGrid, SGrid, GammaGrid);
	
	if(m_log)std::cout<< "Calculate Intensity and autocorrelation for different Energy bins" << std::endl; //std::cin >> dummy;
	/// Calculate Intensity and autocorrelation for different Energy bins
	
	
	TF1* SbdNdS = new TF1((std::string("Integrand S^alpha * dN/dS for ") + source->Name).c_str(),
						[&dNdSSplines] (double* args, double* params) // args[0]: log(S)   params[0]: Gamma  params[1]: alpha  params[2]: intensitybin index
						{ 	//std::cout << "in S^alpha * dN/dS  S = " << args[0] << '\t' << powf(args[0], params[1]) <<  '\t' << source->dNoverdS( args[0], params[0], NULL) <<  std::endl;
							//if(powf(args[0], params[1]
							return exp(args[0]*(1.+params[1])) * std::max(0., dNdSSplines.at((int)params[2])->Eval( exp(args[0]), params[0])); },
						log(SBounds_global.first), log(SBounds_global.second), 3);  			// check range
	
	TF1* GammaIntegrand = new TF1((std::string("Integrated S^alpha * dN/dS for ") + source->Name).c_str(),
										[ SbdNdS, S_tIntensitySpline, &SGrid] (double* args, double* params) // args[0]: Gamma   params[0]: EBin 
										{	SbdNdS->SetParameter(0, args[0]);
											//std::cout << "in Integrated S^alpha * dN/dS: Gamma: " << args[0] << "  S_t = " << S_tSpline->Eval(params[0], args[0]) << std::endl;
											if(S_tIntensitySpline->Eval(params[0], args[0]) <= SGrid.at(0)) return 0.;
											return SbdNdS->Integral(log(SGrid.at(0)), log(S_tIntensitySpline->Eval(params[0], args[0])), 1e-4); },
											source->GammaBounds.first, source->GammaBounds.second, 1);
	
	// first do Intensity
	SbdNdS->SetParameter(1, 1.);  // alpha =1
	
	for(unsigned int i = 0; i < IntensityBins.size(); i++)
	{
		SbdNdS->SetParameter(2, i);
		GammaIntegrand->SetParameter(0, (IntensityBins[i].first + IntensityBins[i].second)/2.);
		
		if(GammaGrid.size() == 2) // No Gamma dependency
		{
			SbdNdS->SetParameter(0, GammaGrid.at(0));
			const double S = S_tIntensitySpline->Eval((IntensityBins.at(i).first + IntensityBins.at(i).second)/2., GammaGrid.at(0));
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
		GammaIntegrand->SetParameter(0, (APSBins.at(i).first + APSBins.at(i).second)/2.);
		
		double* C_p = new double[DT->S_t_1.size()];  // Variate S_t to test dependence on flux threshold
		for(unsigned int j = 0; j < DT->S_t_1.size(); j++)
		{
			if(GammaGrid.size() == 2)   
			{
				const double S = DT->S_t_1.at(j)/DT->S_t_1GeV*S_tAPSSpline->Eval((APSBins.at(i).first + APSBins.at(i).second)/2., GammaGrid.at(0));
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
	// now calc cross correlations
	for(unsigned int i = 0; i < APSBins.size(); i++)
	{
		for(unsigned int j = 0; j < i; j++)
		{
			source->APS->at(i, j) = std::make_shared<gsl1DInterpolationWrapper>(sqrt(*source->APS->at(i, i).get() *  *source->APS->at(j, j).get()));
			source->APS->at(j, i) = source->APS->at(i, j); 
		}
	}*/
	/*
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
	}*/
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


/*
/// Calculates dN/dS 
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtaindNoverdS(AstrophysicalSource* source, Bounds EBin, const std::vector<double>& zGrid, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid)
{
	// Let's us convert S to L in the energy bin
	auto SoverLSpline = ObtainSoverLMapping(source, EBin, zGrid, GammaGrid);
	
	// first calculate inner integral on grid in flux, redshift and photon index
	auto IntegratedrohSpline = std::make_shared<Interpolation3DWrapper>(SGrid.data(), SGrid.size(), zGrid.data(), zGrid.size(), GammaGrid.data(), GammaGrid.size());	
	
	auto Integrandroh = new TF1(("Integrand roh*L/S for" + source->Name).c_str(),
							[SoverLSpline, source] (double* args, double* params) 	// args[0]: S  params[0]: z  params[1]: Gamma
							{	const double SoverL = SoverLSpline->Eval(params[0], params[1]); assert(SoverL<=0); const double L = args[0]/SoverL;
								if(L > source->LumBounds.second || L < source->LumBounds.first) return 0.;
								return source->RescaledLuminosityFunction(L, params[0], params[1])/(SoverL/9.521e48*pow(1._Mpc,2.); },
								SBounds_global.first, SBounds_global.second, 2);
	for(unsigned int k = 0; k < GammaGrid.size(); k++)
	{
		for(unsigned int j = 0; j < zGrid.size(); j++)
		{
			Integrandroh->SetParameters(zGrid.at(j), GammaGrid.at(k));
			for(unsigned int i = 0; i< SGrid.size(); i++)
			{
				IntegratedrohSpline->Val(i, j, k) = Integrandroh->Integral(SGrid.at(i), SGrid.at(i)*1.1, 1e-3);	// cm^2 to Mpc^2
			}
		}
	}
	delete Integrandroh;
	//if(m_log)IntegratedrohSpline->print();
	// now calculate outer integral in redshift
	auto dNdSSpline = std::make_shared<gsl2DInterpolationWrapper>(SGrid.data(), SGrid.size(), GammaGrid.data(), GammaGrid.size());
	
	Bounds zBounds_local; zBounds_local.first = std::max(source->zBounds.first, zBounds_global.first);
	zBounds_local.second = std::min(source->zBounds.second, zBounds_global.second);
	
	auto IntegratedrohdVdz = new TF1(("Integrand dz\\int roh*dV/dz for " + source->Name).c_str(),
										[IntegratedrohSpline, this] (double* args, double* params)	// args[0]:log(z) params[0]: S  params[1]: Gamma
										{ return exp(args[0])*IntegratedrohSpline->Eval(params[0], exp(args[0]), params[1])*CM->ComovingVolumeElement(exp(args[0])); },
										log(zBounds_local.first), log(zBounds_local.second), 2);
	
	for(unsigned int i = 0; i < SGrid.size(); i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			IntegratedrohdVdz->SetParameters(SGrid.at(i), GammaGrid.at(j));
			dNdSSpline->Val(i, j) = IntegratedrohdVdz->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-3)/(0.1*SGrid.at(i));
		}
	}
	delete IntegratedrohdVdz;
	dNdSSpline->Initialize();	

	return dNdSSpline;
}		*/								

/*
 * 
/// Calculates dN/dS 
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtaindNoverdS(AstrophysicalSource* source, Bounds EBin, const std::vector<double>& zGrid, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid)
{
	auto SoverLSpline = ObtainSoverLMapping(source, EBin, zGrid, GammaGrid);
	
	auto dNdSSpline = std::make_shared<gsl2DInterpolationWrapper>(SGrid.data(), SGrid.size(), GammaGrid.data(), GammaGrid.size());
	
	Bounds zBounds_local; zBounds_local.first = std::max(source->zBounds.first, zBounds_global.first);
	zBounds_local.second = std::min(source->zBounds.second, zBounds_global.second);
	
	
	// Since the integral borders of the inner integral depend on the outer integral, we have to do it like this
	
	TF1* rohdVdz = new TF1((std::string("Integrand rho dV/dz for ") + source->Name).c_str(),   
							[source, this] (double *args, double* params) // args[0]: log(L)  params[0]:z   params[1]: Gamma   
							{ 	//std::cout << source->RescaledLuminosityFunction(args[0], params[0], params[1]) << '\t' << CM->ComovingVolume(params[0]) << std::endl;
								//if(isnan(source->RescaledLuminosityFunction(args[0], params[0], params[1]))) std::cout << "nan for L = " << args[0] << "  z = " << params[0] << "  Gamma = " << params[1] << std::endl;
								return  exp(args[0])*source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1])*CM->ComovingVolumeElement(params[0]); },
							log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), 2);
							
	TF1* IntegratedrohdVdzOverL = new TF1((std::string("Integrated(rho dV/dz) over L for ") + source->Name).c_str(),   
										[source, rohdVdz, SoverLSpline] (double *args, double* params) // args[0]: log(z)  params[0]: S  params[1]: Gamma 
										{ 	double SoverL = SoverLSpline->Eval(exp(args[0]), params[1]);
											rohdVdz->SetParameters(exp(args[0]), params[1]);
											std::cout << "S: " << params[0] <<"SOverL: " << SoverL << "   Int bounds: " << params[0] / SoverL << " - " << 1.1*params[0]/SoverL << std::endl;
											return  exp(args[0])*rohdVdz->Integral(log(params[0] / SoverL) , log(1.1*params[0]/SoverL), 1e-4)/1._ergpers; }, // Integrate from log(L(S)) to log(L(S+deltaS)) 
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
}					*/

/// Obtain mapping between flux S and Luminosity L for different redshifts and - if applicable - different photon indeces
std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainSoverLMapping(AstrophysicalSource* source, Bounds EnergyBin, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid)
{
	// InterpolationWrapper that will hold the values				
	auto SoverLSpline = std::make_shared<gsl2DInterpolationWrapper>(zGrid.data(), zGrid.size(), GammaGrid.data(), GammaGrid.size()); 
	
	// Define integrands using ROOT and Lambda functions
	TF1* SIntegrand = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(), 
								[source] (double *args, double *params)   //  args[0]: Energy    params[0]: z   params[1]: Gamma
								{ return source->EnergySpectrum(args[0], params[0], params[1]); },
								/*Emin*/EBounds_global.first, /*Emax*/EBounds_global.second, /*npar*/ 2);
					
	TF1* LIntegrand = new TF1((std::string("Integrand E/k dN/dE for ") + source->Name).c_str(),
								[source] (double *args, double *params)   //  args[0]: Energy    params[0]: z   params[1]: Gamma
								{ return args[0]*source->EnergySpectrumOverK(args[0], params[0], params[1]) ; },
								/*Emin*/   0.1_GeV, /*Emax*/ 100._GeV, /*npar*/ 2);
	
	
	// Integrate
	double L=0, S = 0;
	for(unsigned int i = 0; i < zGrid.size(); i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			SIntegrand->SetParameters(zGrid.at(i), GammaGrid.at(j));
			LIntegrand->SetParameters(zGrid.at(i), GammaGrid.at(j));
			L = LIntegrand->Integral(0.1_GeV , 100._GeV , 1e-3)*1.602e-3  *1._ergpers;		// convert from GeV to erg
			S = SIntegrand->Integral(EnergyBin.first, EnergyBin.second, 1e-3)*1._photonspercm2s;	
			SoverLSpline->Val(i, j) =  S / (L*4.*M_PI* pow(CM->ComovingDistance(zGrid[i]),2.)*9.521e48/pow(1._Mpc,2.)) ;	// Mpc^2 to cm^2
		}
	}
	delete SIntegrand; delete LIntegrand;  				// Clean up	
	
	SoverLSpline->Initialize();// Interpolate S/L over the Grid

	return SoverLSpline;
}

std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > Benchmark::ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<Bounds>& EBins, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline)
{	
	// Calculate flux threshold	
	double** S_t = new double*[EBins.size()]; for(unsigned int i = 0; i < EBins.size(); i++) S_t[i] = new double[GammaGrid.size()];
	
	TF1* dNdE = new TF1((std::string("Integrand dN/dE for ") + source->Name).c_str(),
							[dNdESpline, source] (double *args, double* params) // args[0]: Energy   params[0]: Gamma
							{ return  dNdESpline->Eval(args[0], params[0]); },
							EBounds_global.first, EBounds_global.second, 1);
	
	double denominator = 0;
	for(unsigned int i = 0; i < GammaGrid.size(); i++)
	{
		dNdE->SetParameter(0, GammaGrid.at(i));
		denominator = dNdE->Integral(1._GeV, 100._GeV, 1e-4); 	// We can save some computation time by doing the calculations this way, but need to be careful with pointers
		for(unsigned int j = 0; j < EBins.size(); j++)
		{
			if(denominator <= 0) 
				S_t[j][i] = 0;
			else 	
				S_t[j][i] = dNdE->Integral(EBins.at(j).first, EBins.at(j).second, 1e-4) * DT->S_t_1GeV / denominator;
		} 			// Integrate in specific Energy Bin
		
	}
	delete dNdE;
	std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > S_tSplines;
	for(unsigned int i = 0; i < EBins.size(); i++) 
	{
		S_tSplines.push_back(std::make_shared<gsl1DInterpolationWrapper>(GammaGrid.data(), GammaGrid.size(), S_t[i]));
		delete []S_t[i];
	}
	delete []S_t;

	return S_tSplines;
}

std::shared_ptr<gsl2DInterpolationWrapper> Benchmark::ObtainEffectiveEnergySpectrum(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid)
{
	std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline = ObtainSoverLMapping(source, Bounds(0.1_GeV, 100._GeV), zGrid, GammaGrid);
	
	std::vector<double> dIoverdz; 
	dIoverdz.resize(zGrid.size());
	if(m_log)std::cout<< "Obtaining dI/dz" << std::endl;
	Bounds LBounds_local = Bounds(std::max(source->LumBounds.first, LuminosityBounds_global.first), std::min(source->LumBounds.second, LuminosityBounds_global.second));
	
	/// Obtain dI/dz
	if(GammaGrid.size() == 2)	// no Gamma Dependency
	{
		TF1* Integrand = new TF1("1D Integrand for dI/dz",   
								[source, SoverLSpline, this] (double* args, double* params) // args[0]: log(L)   params[0]: z   params[1]: Gamma, but should be irrelevant
								{ const double S = exp(args[0])* SoverLSpline->Eval(params[0], params[1]);
									//std::cout << "1D Integrand: " << exp(args[0]) << '\t' << S << '\t' << CM->ComovingVolumeElement(params[0]) << '\t' << source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1]) << '\t' << (1.-DT->DetectionEfficiency(S)) << std::endl;
								  return exp(args[0])*S*CM->ComovingVolumeElement(params[0])*source->RescaledLuminosityFunction(exp(args[0]), params[0], params[1])*(1.-DT->DetectionEfficiency(S)); },
								 log(LuminosityBounds_global.first), log(LuminosityBounds_global.second), /*npar*/ 2);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			Integrand->SetParameters(zGrid.at(i), GammaGrid.at(0));
			dIoverdz.at(i) = Integrand->Integral(log(LBounds_local.first), log(LBounds_local.second), 1e-4)/1e48;	// normalization to compare with other
		}		
		delete Integrand;												
	}
	else
	{
		TF2* Integrand = new TF2("2D Integrand for dI/dz",     
								[source, SoverLSpline, this] (double* args, double* params) // args[0]: log(L)  args[1]: Gamma   params[0]: z 
								{ const double S = exp(args[0])* SoverLSpline->Eval(params[0], args[1]);
									return  exp(args[0]) * S * CM->ComovingVolumeElement(params[0]) * source->RescaledLuminosityFunction(exp(args[0]), params[0], args[1]) * (1.-DT->DetectionEfficiency(S)); },
								 log(LuminosityBounds_global.first), log(LuminosityBounds_global.second),
									source->GammaBounds.first, source->GammaBounds.second, /*npar*/ 1);
		for(unsigned int i = 0; i < zGrid.size(); i++)
		{
			Integrand->SetParameters(zGrid.at(i), 0.);
			dIoverdz.at(i) = Integrand->Integral(log(LBounds_local.first), log(LBounds_local.second), 
												source->GammaBounds.first, source->GammaBounds.second , 1e-4)/1e48;
		}		
		delete Integrand;
	}
	
	if(m_log)
	{
		std::cout << "dI/dz = ";
		for(unsigned int i = 0; i < zGrid.size(); i++) std::cout << "(" << zGrid[i] << "," << dIoverdz[i] << ")" << '\t';
	}
	
	/// Use dI/dz to calculate effective energy spectrum
	auto dIdzSpline = std::make_shared<gsl1DInterpolationWrapper>(zGrid.data(), zGrid.size(), dIoverdz.data(), gsl_interp_linear, 0); 
	
	if(m_plot)
	{
		auto g = dIdzSpline->MakeGraph();  g->SetName(source->Name.c_str());
		dIdzFile->cd();
		g->Write();
	}	
	
	if(m_log) std::cout<< std::endl << "calculate effective energy spectrum" << std::endl; //std::cin >> dummy;
	Bounds zBounds_local; zBounds_local.first = zGrid[0]; zBounds_local.second = zGrid[zGrid.size()-1];
	
	auto dNdESpline = std::make_shared<gsl2DInterpolationWrapper>(EGrid.data(), EGrid.size(), GammaGrid.data(), GammaGrid.size());
	
	// Integrate it to get dN/dE without z dependence
	
	TF1* dIdzdNdE = new TF1((std::string("Integrand dI/dz * dN/dE for ") + source->Name).c_str(),
							[dIdzSpline, source] (double* args, double* params) // args[0]: log z  params[0]:E  params[1]:Gamma
							{ //std::cout << "z = " << exp(args[0]) << "  E = " << params[0] << "  Gamma = " << params[1]
								//	<< "   dI/dz = " << dIdzSpline->Eval(exp(args[0])) << "  dN/dE = " << source->EnergySpectrum(params[0], args[0], params[1]) << std::endl;
							 return exp(args[0])* dIdzSpline->Eval(exp(args[0])) * source->EnergySpectrum(params[0], exp(args[0]), params[1]); },
							log(zBounds_local.first), log(zBounds_local.second), /*npar*/ 2);
							
							
	TF1* dIdz = new TF1((std::string("Integrand dI/dz for ") + source->Name).c_str() ,
							[dIdzSpline, source] (double* args, double* params) // args[0]: log(z) 
							{ return exp(args[0])*dIdzSpline->Eval(exp(args[0])); },
							log(zBounds_local.first), log(zBounds_local.second), /*npar*/ 0);
							
	double denominator = dIdz->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-4);
	for(unsigned int i = 0; i < EGrid.size() ; i++)
	{
		for(unsigned int j = 0; j < GammaGrid.size(); j++)
		{
			if(denominator <= 0) dNdESpline->Val(i, j) = 0;
			else
			{
				dIdzdNdE->SetParameters(EGrid.at(i), GammaGrid.at(j));
				dNdESpline->Val(i, j) = dIdzdNdE->Integral(log(zBounds_local.first), log(zBounds_local.second), 1e-3)/denominator;
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
	assert(zGridLen >=2); assert(kGridLen >=2);
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
	auto P1HaloIntegrand = new TF1(("P1 Halo Integrand for " + DM->Name).c_str(),
									[this, DM] (double* args, double* params)	// args[0]: log(M)  params[0]: k   params[1]: z
									{	return exp(args[0])* HM->HaloMassFunction(exp(args[0]), params[1]) * pow(DM->SourceDensityFT(params[0], exp(args[0]), params[1]), 2.);	},
									log(HM->MBounds.first), log(HM->MBounds.second), 2);
	
	auto P2HaloIntegrand = new TF1(("The Integrand dn/dm * LinearHaloBias *sourcedensityFT for the 2 Halo term for " + DM->Name).c_str(),
									[DM, this] (double* args, double* params) // args[0]: log(M)  params[0]: k  params[1]: z
									{ //std::cout << "P2: " << exp(args[0]) << '\t' << params[0] << '\t' << HM->HaloMassFunction(exp(args[0]), params[1]) << '\t' << HM->LinearHaloBias(exp(args[0]), params[1]) << '\t' << DM->SourceDensityFT(params[0], exp(args[0]), params[1]) << std::endl;
										return exp(args[0])* HM->HaloMassFunction(exp(args[0]), params[1]) * HM->LinearHaloBias(exp(args[0]), params[1]) * DM->SourceDensityFT(params[0], exp(args[0]), params[1]); },
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
			
			//if(m_log) std::cout << "k = " << kGrid.at(i) << "  z = " << zGrid.at(j) << "  p1 = " << p1 << "  p2 = " << p2 << std::endl;
			//if(m_log) std::cout << '(' << i << ", " << j << "): " << _3DPowerSpectrumSpline->Val(i, j) << std::endl;
		}
	} 
	delete P1HaloIntegrand; delete P2HaloIntegrand;
	_3DPowerSpectrumSpline->Initialize();
	if(m_log) {std::cout << "3D Power Spectrum: "; _3DPowerSpectrumSpline->print(); }
	
	if(m_plot)
	{
		std::vector<int> z = {0, 1, 3};
		_3DPSFile->cd();
		for(unsigned int i = 0; i < z.size(); i++)
		{
			auto g = _3DPowerSpectrumSpline->MakeGraphAlongX(z.at(i)); g->SetName((DM->Name + std::to_string(z.at(i))).c_str());		
			g->Write();
		}
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
				DM->APS->at(i, j, k) = APSCrossIntegrand->Integral(log(zBounds_global.first), log(zBounds_global.second), EBins.at(i).first, EBins.at(i).second,
															EBins.at(j).first, EBins.at(j).second, 1e-3);
				DM->APS->at(j, i, k) = DM->APS->at(i, j, k);
				//if(m_log) std::cout << "( " << i << ", " << j << ", " << k << "): " << DM->APS->at(i,j,k) << std::endl;
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
			APSAutoIntegrand->SetParameter(0, Multipoles.at(k));
			DM->APS->at(i, i, k) = APSAutoIntegrand->Integral(log(zBounds_global.first), log(zBounds_global.second), EBins.at(i).first, EBins.at(i).second, 1e-3);
			//if(m_log) std::cout << "( " << i << ", " << k << "): " << DM->APS->at(i,i,k) << std::endl;
		}
	}
	delete APSAutoIntegrand;
}

#endif
