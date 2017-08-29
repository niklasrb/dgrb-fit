#ifndef DGRBFIT_H
#define DGRBFIT_H

#include "multinest.h"
#include "Constants.h"
#include "AstrophysicalSource.h"
#include "DarkMatter.h"

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
				APS->at(i,j) = std::make_shared<gsl2DInterpolationWrapper>(gsl2DInterpolationWrapper::Combine(APSdata, parameter));		// first parameter is S_t, second another
			}
		}
	}
	
	double ScaleParameter(double& param)
	{
		return parameterBounds.first + (parameterBounds.second - parameterBounds.first)*param;
	}
};





class IntensityFit : public MultiNestFit
{
protected:
	std::vector<Bounds> IntensityBins;
	std::vector<Measurement> IntensityMeasurement;
	
	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::vector<std::shared_ptr<DarkMatter> > dmModels;
	std::vector<std::shared_ptr<AstrophysicalSourceClass> > AstrophysicalSourceClasses;
	
	
	double IntensityChiSquared(double* parameters)
	{
		double chiSquared = 0, intensity;
		for(unsigned int i = 0; i < IntensityBins.size(); i++)
		{
			intensity = 0;
			for(unsigned int j =0; j < AstrophysicalSources.size(); j++)
			{
				intensity += pow(10, 2*parameters[j] - 1)*AstrophysicalSources[j]->Intensity[i];	// parameter is coefficient
			}
			int n = AstrophysicalSources.size();
			for(unsigned int j = 0; j < dmModels.size(); j++)
			{
				intensity += pow(10, 10.*parameters[n+j] - 5.)*dmModels[j]->Intensity[i];
			}
			n += dmModels.size();
			for(unsigned int j = 0; j < AstrophysicalSourceClasses.size(); j++)
			{
				auto asc = AstrophysicalSourceClasses[j];
				intensity += pow(10, 2*parameters[n+2*j]-1)*asc->Intensity[i]->Eval(asc->parameterBounds.first + (asc->parameterBounds.second - asc->parameterBounds.first)*parameters[n+2*j+1]);	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
			}
			chiSquared += pow( (IntensityMeasurement.at(i).first - intensity)/IntensityMeasurement.at(i).second, 2);
		}
		return chiSquared/IntensityBins.size();
	}
	
	double loglike(double* Cube, int *npar) override
	{
		return -0.5*IntensityChiSquared(Cube);
	}
	
	void dumper(MultiNestFitData& data) override
	{
		
		std::cout << "Last parameters: "; for(int i = 0; i < data.nPar; i++) std::cout << i << ": "<< data.lastParameters[i] << (i < data.nPar-1 ? '\t' : '\n');
		std::cout << "loglike: " << data.logLike << std::endl;
	}

public:
	IntensityFit(const std::vector<Bounds>& IntensityBins, const std::vector<Measurement>& IntensityMeasurement, const std::vector<std::shared_ptr<AstrophysicalSource> >& AstrophysicalSources, const std::vector<std::shared_ptr<DarkMatter> >& dmModels,
					const std::vector<std::shared_ptr<AstrophysicalSourceClass> >& AstrophysicalSourceClasses, std::string output = "") 
	: MultiNestFit(output) , IntensityBins(IntensityBins), IntensityMeasurement(IntensityMeasurement), AstrophysicalSources(AstrophysicalSources), dmModels(dmModels), AstrophysicalSourceClasses(AstrophysicalSourceClasses)
	{
		assert(IntensityBins.size() == IntensityMeasurement.size() && IntensityBins.size() > 0);
		this->nParameters = AstrophysicalSources.size() + dmModels.size()  + 2* AstrophysicalSourceClasses.size();
	}
	
	void printResults()
	{
		std::cout << "Intensity Fit: " << std::endl;
		for(unsigned int i = 0; i < AstrophysicalSources.size(); i++)
			std::cout << AstrophysicalSources[i]->Name << ": " << pow(10, 2*data.lastParameters[i]-1) << std::endl;
		int n = AstrophysicalSources.size();
		for(unsigned int i = 0; i < dmModels.size(); i++)
			std::cout << dmModels[i]->Name << ": " << pow(10, 10.*data.lastParameters[i+n] -5.) << std::endl;
		for(unsigned int i = 0; i < AstrophysicalSourceClasses.size(); i++)
			std::cout << AstrophysicalSourceClasses[i]->Name << ": " << pow(10, 2*data.lastParameters[n+i]-1) << ",  " << AstrophysicalSourceClasses[i]->ScaleParameter(data.lastParameters[n+i+1]) << std::endl;
	}
	
};


class IntensityAndAPSFit : public IntensityFit
{
protected:
	std::vector<Bounds> APSBins;
	std::shared_ptr<AngularPowerSpectrum<Measurement> > APSMeasurement;
	
	Bounds SBounds;
	
	double APSChiSquared(double* parameters)
	{
		double chiSquared = 0, val;
		//int nS = nPar; 	// position of S parameter
		double S_t = SBounds.first + (SBounds.second - SBounds.first)*parameters[nParameters-1];	// S_t is last parameter
		std::cout << "APSChiSquared:  S = " << S_t <<std::endl;
		for(unsigned int i = 0; i < APSMeasurement->Bin1Size(); i++)
		{
			for(unsigned int j = 0; j <= i; j++)
			{
				for(unsigned int k = 0; k < APSMeasurement->MultipoleNumber(); k++)
				{
					val = 0;
					for(unsigned int l =0; l < AstrophysicalSources.size(); l++)
					{
						val += pow(10, (2.*parameters[l] - 1.))*AstrophysicalSources[l]->APS->at(i, j)->Eval(S_t);	// parameter is coefficient
					}
					int n = AstrophysicalSources.size();
					for(unsigned int l = 0; l < dmModels.size(); l++)
					{
						val += pow(10, 2.*(10.*parameters[n+l] -5.))* dmModels[l]->APS->at(i, j, k);		// this dependence is squared
					}
					n += dmModels.size();
					for(unsigned int l = 0; l < AstrophysicalSourceClasses.size(); l++)
					{
						auto asc = AstrophysicalSourceClasses[l];
						val += pow(10, (2.*parameters[n+2*l]-1.))*asc->APS->at(i, j)->Eval(S_t, asc->ScaleParameter(parameters[n+2*l+1]));	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
					}
					
					chiSquared += pow( (APSMeasurement->at(i, j, k).first - val)/APSMeasurement->at(i, j, k).second, 2);
				}
			}
		}
		return chiSquared/(pow(APSBins.size(),2)*APSMeasurement->MultipoleNumber()/2.);
	}
	
	double loglike(double* Cube, int *npar) override
	{
		return -0.5*(APSChiSquared(Cube) + IntensityChiSquared(Cube));
	}
	
	void dumper(MultiNestFitData& data) override
	{
		
		std::cout << "Last parameters: "; for(int i = 0; i < data.nPar; i++) std::cout << i << ": "<< data.lastParameters[i] << (i < data.nPar-1 ? '\t' : '\n');
		std::cout << "loglike: " << data.logLike << std::endl;
	}
	
public:
	IntensityAndAPSFit(const std::vector<Bounds>& IntensityBins, const std::vector<Measurement>& IntensityMeasurement, const std::vector<Bounds>& APSBins, const std::shared_ptr<AngularPowerSpectrum<Measurement> >& APSMeasurement,
		const std::vector<std::shared_ptr<AstrophysicalSource> >& AstrophysicalSources, const std::vector<std::shared_ptr<DarkMatter> >& dmModels, const std::vector<std::shared_ptr<AstrophysicalSourceClass> >& AstrophysicalSourceClasses, Bounds SBounds, std::string output = "")
		: IntensityFit(IntensityBins, IntensityMeasurement, AstrophysicalSources, dmModels, AstrophysicalSourceClasses, output) , APSBins(APSBins), APSMeasurement(APSMeasurement), SBounds(SBounds)
		{
			assert(APSBins.size() == APSMeasurement->Bin1Size()  && APSBins.size() == APSMeasurement->Bin2Size() && APSBins.size() >=2);
			
			this->nParameters = AstrophysicalSources.size() + dmModels.size() + 2* AstrophysicalSourceClasses.size() + 1; 
			options.nlive = 100;
		} 
	
	void printResults()
	{
		IntensityFit::printResults();
		std::cout << "S_t = " << SBounds.first + (SBounds.second - SBounds.first)*data.lastParameters[nParameters-1] << std::endl;
	}
	
};



#endif
