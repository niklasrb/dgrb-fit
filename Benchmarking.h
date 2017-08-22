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
#include <cstdlib>
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
#include "InterpolationWrapper.h"
#include "EBLAbsorbtionCoefficient.h"
#include "CosmologyModel.h"
#include "HaloModel.h"
#include "Source.h"
#include "AstrophysicalSource.h"
#include "DarkMatter.h"
#include "GalaxyCatalog.h"
#include "AngularPowerSpectrum.h"
#include "Detector.h"
#include "CanvasWrapper.h"


class Benchmark
{
protected:
	bool m_log;
	bool m_plot;
	std::shared_ptr<Canvas> dIdzCanvas;
	//std::shared_ptr<TCanvas> dNdSCanvas;
	
public:
	Bounds LuminosityBounds_global;
	Bounds zBounds_global;  unsigned int zGridLen;
	Bounds kBounds_global;  unsigned int kGridLen;
	Bounds EBounds_global;  unsigned int EGridLen;
	Bounds SBounds_global;  unsigned int SGridLen;
	unsigned int GammaGridLen;
	
	std::vector<Bounds> IntensityBins;
	std::vector<Bounds> APSBins;
	
private:
	std::shared_ptr<CosmologyModel> CM;
	std::shared_ptr<HaloModel> HM;
	std::shared_ptr<Detector> DT;
	//std::vector<Bounds> EBins;
	
public:	
	
	Benchmark(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<Detector> DT, bool log, bool plot) ;
	~Benchmark();
	
	void calculateIntensityAndAPSForAstrophysicalSources(std::vector<std::shared_ptr<AstrophysicalSource> > sources);
	
	void calculateIntensityAndAutocorrelationForDM(std::vector<std::shared_ptr<DarkMatter> > DM);
	
private:
	// For Astrophysical Sources
	void calculateIntensityAndAPS(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& SGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainSoverLMapping(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtaindNoverdS(AstrophysicalSource* source, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainEffectiveEnergySpectrum(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<Bounds>& EGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> SoverLSpline, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline);
	void plotIntermediates(AstrophysicalSource* source, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline);
	
	// For DM
	void calculateAPSForDM(std::shared_ptr<DarkMatter> DM, const std::vector<Bounds>& EBins, const std::vector<double>& zGrid, const std::vector<double>& kGrid, const std::vector<int>& Multipoles);
	void calculateIntensityForDM(std::shared_ptr<DarkMatter> DM, const std::vector<Bounds>& EBins);
};


Benchmark::Benchmark(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<Detector> DT, bool log = true, bool plot = false)  
												: m_log(log), m_plot(plot), CM(CM), HM(HM), DT(DT)
{
	if(m_plot)
	{
		dIdzCanvas = std::make_shared<Canvas>("dIdz", "dI/dz " );
		(*dIdzCanvas)().SetLogx(1); (*dIdzCanvas)().SetLogy(1); 
	//	//dNdSCanvas = std::make_shared<TCanvas>("dNdS", "dN/dS");
	}
}


Benchmark::~Benchmark()
{
}

#include "Benchmarking.cpp"

#endif
