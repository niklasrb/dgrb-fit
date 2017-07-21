#ifndef DGRBFIT_H
#define DGRBFIT_H

#include "multinest.h"
#include "Constants.h"
#include "AstrophysicalSource.h"

#include <tuple>
#include <vector>



#define IntensityMeasurement std::pair< double, double>   // value and error

class IntensityFit : public MultiNestFit
{
	std::vector< Bounds > EBins;
	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::vector< IntensityMeasurement> Intensities;
	
	
	double IntensityChiSquared(double* coefficients)
	{
		double chiSquared = 0;
		for(unsigned int i = 0; i < EBins.size(); i++)
		{
			double intensity = 0;
			for(unsigned int j =0; j < AstrophysicalSources.size(); j++)
			{
				intensity += coefficients[j]*AstrophysicalSources.at(j)->Intensity.at(i).second;
			}
			chiSquared += pow( (Intensities.at(i).first - intensity)/Intensities.at(i).second, 2);
		}
		return chiSquared/EBins.size();
	}
	
	double loglike(double* Cube) override
	{
		return -0.5*IntensityChiSquared(Cube);
	}
	
	IntensityFit(std::vector<Bounds> EBins, std::vector<std::shared_ptr<AstrophysicalSource>> AstrophysicalSources) 
	: MultiNestFit(AstrophysicalSources.size()), EBins(EBins), AstrophysicalSources(AstrophysicalSources) {}
};



#endif
