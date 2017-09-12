#ifndef EBLABSORBTIONCOEFFICIENT_H
#define EBLABSORBTIONCOEFFICIENT_H

#include <cmath>
#include "InterpolationWrapper.h"

#include "TROOT.h"
#include "TFile.h"


/* This class models the EBL absorption coeffiecient
 * I never know how to spell Absor(b/p)tion, sorry for that
 * Just holds a 2D interpolation object
 * But could be changed to whatever model
 */
class EBLAbsorbtionCoefficient
{
protected:
	gsl2DInterpolationWrapper interp;

public:

	EBLAbsorbtionCoefficient(gsl2DInterpolationWrapper interp) : interp(interp) {}

	double operator ()(const double E, const double z)
	{
		return interp.Eval(E, z);
	}
	
	void plot(std::string file)
	{
		auto f = new TFile(file.c_str(), "RECREATE");
		std::vector<double> z = {0.01, 0.1, 2., 10.};
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
