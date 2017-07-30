#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include <vector>
#include <cassert> 
#include <ctime>
#include <ncurses.h>
#include <memory>
#include <limits>
#include <algorithm>
#include <fstream>
#include "TROOT.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/IntegratorOptions.h"

#include "Constants.h"
#include "gsl2DInterpolationWrapper.h"
#include "EBLAbsorbtionCoefficient.h"
#include "CosmologyModel.h"
#include "HaloModel.h"
#include "Source.h"
#include "AstrophysicalSource.h"
#include "DarkMatter.h"
#include "GalaxyCatalog.h"
#include "AngularPowerSpectrum.h"


class Benchmark
{
protected:
	bool m_log;
	//bool m_plot;
	//std::shared_ptr<TCanvas> dIdzCanvas;
	//std::shared_ptr<TCanvas> dNdSCanvas;
	
	Bounds LuminosityBounds_global;
	Bounds zBounds_global;
	Bounds kBounds_global;
	
	std::shared_ptr<CosmologyModel> CM;
	std::shared_ptr<HaloModel> HM;
	std::vector<Bounds> EBins;
	
public:	
	
	Benchmark(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<HaloModel> HM, std::vector<Bounds> EBins, Bounds LBounds, Bounds zBounds, Bounds kBounds, bool log) ;
	~Benchmark();
	
	void calculateIntensityAndAutocorrelationForAstrophysicalSources(std::vector<std::shared_ptr<AstrophysicalSource> > sources, int zGridLen, int GammaGridLen);
	
	void calculateIntensityAndAutocorrelationForDM(std::vector<std::shared_ptr<DarkMatter> > DM, unsigned int zGridLen, unsigned int kGridLen);
	
private:
	// For Astrophysical Sources
	void calculateIntensityAndAutocorrelation(AstrophysicalSource* source,const  std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainSoverLMapping(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	void ObtaindNoverdS(AstrophysicalSource* source, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline, const double S_t_1Gev);
	
	// For DM
	void calculateIntensityAndAutocorrelationForDM(std::shared_ptr<DarkMatter> DM, const std::vector<double>& zGrid, const std::vector<double>& kGrid, const std::vector<int>& Multipoles);
	
};
Benchmark::Benchmark(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<HaloModel> HM, std::vector<Bounds> EBins, Bounds LBounds, Bounds zBounds, Bounds kBounds, bool log = true)  
																			: m_log(log)/*, m_plot(plot)*/,  LuminosityBounds_global(LBounds), zBounds_global(zBounds), kBounds_global(kBounds), CM(CM), HM(HM), EBins(EBins)
{
	//if(m_plot)
	//{
	//	dIdzCanvas = std::make_shared<TCanvas>("dIdz", "dI/dz"); dIdzCanvas->Draw();
	//	//dNdSCanvas = std::make_shared<TCanvas>("dNdS", "dN/dS");
	//}
}


Benchmark::~Benchmark()
{
}

#include "Benchmarking.cpp"

#endif
