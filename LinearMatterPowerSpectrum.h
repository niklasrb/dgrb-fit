#ifndef LINEARMATTERPOWERSPECTRUM_H
#define LINEARMATTERPOWERSPECTRUM_H


#include "InterpolationWrapper.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"

/* Wrapper class that holds the data for the linear matter power spectrum
 * just a 2D interpolation object
 */

class LinearMatterPowerSpectrum
{
protected:
	gsl2DInterpolationWrapper interp;
	
public:
	LinearMatterPowerSpectrum(gsl2DInterpolationWrapper interp) : interp(interp) {}
	
	double operator ()(const double k, const double z)
	{
		return interp.Eval(k*1._Mpc, z)/pow(1._Mpc,3);
	}
	
	void plot(std::string file)
	{
		auto f = new TFile(file.c_str(), "RECREATE");
		std::vector<double> z = {1., 5., 10., 20.};
		for(unsigned int i = 0; i < z.size(); i++)
		{
			auto g = interp.MakeGraphAlongX(z[i]);
			g->SetName(std::to_string(z[i]).c_str());
			g->Write();
		}
		f->Close();
	}
};


#endif
