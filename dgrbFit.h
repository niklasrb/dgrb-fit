#ifndef DGRBFIT_H
#define DGRBFIT_H

#include "multinest.h"
#include "Constants.h"
#include "AstrophysicalSource.h"

#include <tuple>
#include <vector>
#include <memory>


struct AstrophysicalSourceClass
{
	std::string Name;
	std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > Intensity;		// If the intensities depend on another paramter, i.e. E_cut
	std::shared_ptr<AstrophysicalSourceAPS<std::shared_ptr<gsl2DInterpolationWrapper> > > APS;	// Depending on S_t_1Gev and a second parameter 
	Bounds parameterBounds;
	
	AstrophysicalSourceClass(const std::vector<Bounds>& IntensityBins, const std::vector<Bounds>& APSBins, const std::vector<std::shared_ptr<AstrophysicalSource> >& sources, const std::vector<double>& parameter)
	{
		assert(sources.size() >= 2); 
		assert(sources.size() == parameter.size());
		// Merge intensity data
		std::vector<double> Intensities; Intensities.resize(parameter.size());		
		for(unsigned int i = 0; i < IntensityBins.size(); i++)
		{
			for(unsigned int j = 0; j < sources.size(); j++) Intensities[j] = sources[j]->Intensity[i];
			Intensity.push_back(std::make_shared<gsl1DInterpolationWrapper>(parameter.data(), parameter.size(), Intensities.data()));
		}
		// Merge APS data
		std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > APSdata; APSdata.resize(parameter.size());
		for(unsigned int i = 0; i < APSBins.size(); i++)
		{
			for(unsigned int j = 0; j < APSBins.size(); j++)
			{
				for(unsigned int k = 0; k < parameter.size(); k++) APSdata[k] = sources[k]->APS->at(i,j);
				APS->at(i,j) = std::make_shared<gsl2DInterpolationWrapper>(gsl2DInterpolationWrapper::Combine(APSdata, parameter));
			}
		}
	}
};


typedef  std::pair< double, double> IntensityMeasurement;  // value and error


class IntensityFit : public MultiNestFit
{
protected:
	std::vector< Bounds > EBins;
	std::vector< IntensityMeasurement> Measurement;
	std::vector<std::shared_ptr<DGRBSource> > Sources;
	std::vector<std::shared_ptr<AstrophysicalSourceClass> > AstrophysicalSourceClasses;
	
	
	double IntensityChiSquared(double* parameters)
	{
		double chiSquared = 0;
		for(unsigned int i = 0; i < EBins.size(); i++)
		{
			double intensity = 0;
			for(unsigned int j =0; j < Sources.size(); j++)
			{
				intensity += parameters[j]*Sources[j]->Intensity[i];	// parameter is coefficient
			}
			int n = Sources.size();
			for(unsigned int j = 0; j < AstrophysicalSourceClasses.size(); j++)
			{
				intensity += parameters[n+j]*AstrophysicalSourceClasses[j]->Intensity[i]->Eval(parameters[n+j+1]);	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
			}
			chiSquared += pow( (Measurement.at(i).first - intensity)/Measurement.at(i).second, 2);
		}
		return chiSquared/EBins.size();
	}
	
	double loglike(double* Cube, int *npar) override
	{
		return -0.5*IntensityChiSquared(Cube);
	}

public:
	IntensityFit(std::vector<Bounds> EBins, std::vector<IntensityMeasurement> Measurement, const std::vector<std::shared_ptr<DGRBSource> >&Sources, std::vector<std::shared_ptr<AstrophysicalSourceClass> > AstrophysicalSourceClasses) 
	: MultiNestFit() , EBins(EBins), Measurement(Measurement), Sources(Sources), AstrophysicalSourceClasses(AstrophysicalSourceClasses)
	{
		this->nParameters = Sources.size() + 2* AstrophysicalSourceClasses.size();
	}
	
};



#endif
