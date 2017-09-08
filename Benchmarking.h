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
#include "TFile.h"
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

	std::shared_ptr<TFile> dIdzFile;
	std::shared_ptr<TFile> NofSFile;
	std::shared_ptr<TFile> _3DPSFile;
	std::shared_ptr<TFile> IntensityFile;
	std::shared_ptr<TFile>  dNdEFile;
	
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
	
public:	
	
	Benchmark(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<Detector> DT, bool log, bool plot, std::string) ;
	~Benchmark();
	
	void calculateIntensityAndAPSForAstrophysicalSources(std::vector<std::shared_ptr<AstrophysicalSource> > sources);
	void calculateIntensityAndAutocorrelationForDM(std::vector<std::shared_ptr<DarkMatter> >& DM, const std::vector<double>& Multipoles);
	
private:
	// For Astrophysical Sources
	void calculateIntensityAndAPS(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& SGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainSoverLMapping(AstrophysicalSource* source, Bounds EnergyBin, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtaindNoverdS(AstrophysicalSource* source, Bounds EnergyBin, const std::vector<double>& zGrid, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid);
	std::shared_ptr<gsl2DInterpolationWrapper> ObtainEffectiveEnergySpectrum(AstrophysicalSource* source, const std::vector<double>& EGrid, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<Bounds>& EBins, const std::vector<double>& GammaGrid, std::shared_ptr<gsl2DInterpolationWrapper> dNdESpline);
	
	// For DM
	void calculateAPSForDM(std::shared_ptr<DarkMatter>& DM, const std::vector<Bounds>& EBins, const std::vector<double>& zGrid, const std::vector<double>& kGrid, const std::vector<double>& Multipoles);
	void calculateIntensityForDM(std::shared_ptr<DarkMatter> DM, const std::vector<Bounds>& EBins);
};


Benchmark::Benchmark(std::shared_ptr<CosmologyModel> CM, std::shared_ptr<HaloModel> HM, std::shared_ptr<Detector> DT, bool log = true, bool plot = false, std::string dir = "")  
												: m_log(log), m_plot(plot), CM(CM), HM(HM), DT(DT)
{
	if(m_plot)
	{
		dIdzFile = std::make_shared<TFile>((dir + "dIdz.root").c_str(), "RECREATE");
		NofSFile = std::make_shared<TFile>((dir + "NofS.root").c_str(), "RECREATE");
		_3DPSFile = std::make_shared<TFile>((dir + "3DPS.root").c_str(), "RECREATE");
		IntensityFile = std::make_shared<TFile>((dir + "Intensity.root").c_str(), "RECREATE");
		dNdEFile = std::make_shared<TFile>((dir + "dNdE.root").c_str(), "RECREATE");
	}
}



#include "Benchmarking.cpp"

#endif
