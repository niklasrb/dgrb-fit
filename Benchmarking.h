#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include <vector>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <cassert> 
#include "TROOT.h"
#include "TF1.h"
#include "TF2.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/IntegratorOptions.h"
#include "Source.h"
#include "AstrophysicalSource.h"
#include "Constants.h"
#include "CosmologyModel.h"
#include "Detector.h"
#include <ctime>
#include <ncurses.h>
#include <memory>
#include <limits>
#include <algorithm>

class Benchmark
{
private:
	bool m_log;
	bool m_plot;
	std::shared_ptr<TCanvas> dIdzCanvas;
	std::shared_ptr<TCanvas> dNdSCanvas;
	Bounds LuminosityBounds;
	
	
public:
	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::shared_ptr<CosmologyModel> CM;
	std::shared_ptr<Detector> D;
	std::vector<Bin> EBins;
	
	Benchmark(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<Detector> _D, bool log, bool plot);
	~Benchmark();
	
	void calculateIntensityAndAutocorrelation(int zGridLen, int GammaGridLen);
	
private:
	// For Astrophysical Sources
	void calculateIntensityAndAutocorrelation(AstrophysicalSource* source,const  std::vector<double>& zGrid, const std::vector<double>& GammaGrid);
	void ObtainSoverLMapping(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl_spline2d>& SoverLSpline, std::shared_ptr<gsl_interp_accel>& zAcc, std::shared_ptr<gsl_interp_accel>& GammaAcc);
	void ObtaindNoverdS(AstrophysicalSource* source, const std::vector<double>& SGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl_spline2d>& SoverLSpline, std::shared_ptr<gsl_interp_accel>& zAcc, std::shared_ptr<gsl_interp_accel>& GammaAcc);
	std::vector<double> ObtainFluxThreshold(AstrophysicalSource* source, const std::vector<double>& zGrid, const std::vector<double>& GammaGrid, std::shared_ptr<gsl_spline2d>& SoverLSpline, std::shared_ptr<gsl_interp_accel>& zAcc, std::shared_ptr<gsl_interp_accel>& GammaAcc);
};
Benchmark::Benchmark(std::shared_ptr<CosmologyModel> _CM, std::shared_ptr<Detector> _D, bool log, bool plot )  : m_log(log), m_plot(plot),  CM(_CM), D(_D)
{
	LuminosityBounds.first = 1e1; LuminosityBounds.second = 1e15;
	if(m_plot)
	{
		dIdzCanvas = std::make_shared<TCanvas>("dIdz", "dI/dz");
		dNdSCanvas = std::make_shared<TCanvas>("dNdS", "dN/dS");
	}
}


Benchmark::~Benchmark()
{
}

#include "Benchmarking.cpp"

#endif
