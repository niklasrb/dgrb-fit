#ifndef DGRBFIT_H
#define DGRBFIT_H

#include "multinest.h"
#include "Constants.h"
#include "AstrophysicalSource.h"
#include "DarkMatter.h"

#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include <tuple>
#include <vector>
#include <memory>


struct AstrophysicalSourceClass
{
	std::string Name;
	std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > Intensity;		// If the intensities depend on another paramter, i.e. E_cut
	std::shared_ptr<AstrophysicalSourceAPS<std::shared_ptr<gsl2DInterpolationWrapper> > > APS;	// Depending on S_t_1Gev and a second parameter 
	Bounds parameterBounds;
	
	AstrophysicalSourceClass(const std::vector<Bounds>& IntensityBins, const std::vector<Bounds>& APSBins, const std::vector<std::shared_ptr<AstrophysicalSource> >& sources, const std::vector<double>& parameter, std::string name = "")
	{
		assert(sources.size() >= 2); 
		assert(sources.size() == parameter.size());
		if(name.empty()) Name = sources.at(0)->Name;
		else Name = name;
		// Take Bounds
		parameterBounds = Bounds(parameter.at(0), parameter.at(parameter.size()-1));
		// Merge intensity data
		std::vector<double> Intensities; Intensities.resize(parameter.size());		
		for(unsigned int i = 0; i < IntensityBins.size(); i++)
		{
			for(unsigned int j = 0; j < sources.size(); j++) Intensities[j] = sources[j]->Intensity[i];
			Intensity.push_back(std::make_shared<gsl1DInterpolationWrapper>(parameter.data(), parameter.size(), Intensities.data()));
		}
		// Merge APS data
		APS = std::make_shared<AstrophysicalSourceAPS<std::shared_ptr<gsl2DInterpolationWrapper> > >(APSBins);
		std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > APSdata; APSdata.resize(parameter.size());
		for(unsigned int i = 0; i < APSBins.size(); i++)
		{
			for(unsigned int j = 0; j < APSBins.size(); j++)
			{
				for(unsigned int k = 0; k < parameter.size(); k++) APSdata.at(k) = sources.at(k)->APS->at(i,j);
				APS->at(i,j) = std::make_shared<gsl2DInterpolationWrapper>(gsl2DInterpolationWrapper::Combine(APSdata, parameter));		// first parameter is S_t, second another
			}
		}
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
	
	
	double IntensityChiSquared()
	{
		double chiSquared = 0, intensity;
		for(unsigned int i = 0; i < IntensityBins.size(); i++)
		{
			intensity = 0;
			for(unsigned int j =0; j < AstrophysicalSources.size(); j++)
			{
				intensity += pow(10, Parameters.at(j)())*AstrophysicalSources.at(j)->Intensity.at(i);	// parameter is coefficient
			}
			int n = AstrophysicalSources.size();
			for(unsigned int j = 0; j < dmModels.size(); j++)
			{
				intensity += pow(10, Parameters.at(n+j)())*dmModels.at(j)->Intensity.at(i);
			}
			n += dmModels.size();
			for(unsigned int j = 0; j < AstrophysicalSourceClasses.size(); j++)
			{
				auto asc = AstrophysicalSourceClasses.at(j);
				intensity += pow(10, Parameters.at(n+2*j)())*asc->Intensity.at(i)->Eval(Parameters.at(n+2*j+1)());	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
			}
			chiSquared += pow( (IntensityMeasurement.at(i).first - intensity)/IntensityMeasurement.at(i).second, 2);
		}
		return chiSquared/IntensityBins.size();
	}
	
	double loglike() override
	{
		return -0.5*IntensityChiSquared();
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
		for(unsigned int i = 0; i < AstrophysicalSources.size(); i++) Parameters.push_back(MultiNestParameter(Bounds(3, 6), AstrophysicalSources[i]->Name + " LumFunc log factor"));
		for(unsigned int i = 0; i < dmModels.size(); i++) Parameters.push_back(MultiNestParameter(Bounds(-5, 5), dmModels[i]->Name + " Windowfunction log factor"));
		for(unsigned int i = 0; i < AstrophysicalSourceClasses.size(); i++)
		{
			Parameters.push_back(MultiNestParameter(Bounds(3, 6), AstrophysicalSourceClasses[i]->Name + " LumFunc log factor"));
			Parameters.push_back(MultiNestParameter(AstrophysicalSourceClasses[i]->parameterBounds, AstrophysicalSourceClasses[i]->Name + " second param"));
		}
	}
	
	void printResults()
	{
		std::cout << "Intensity Fit: " << std::endl;
		for(unsigned int i = 0; i < Parameters.size(); i++) Parameters.at(i).print();
		std::cout << "chi^2 = " << IntensityChiSquared() << std::endl;
	}
	
	void plotResults(TFile* f)
	{
		f->cd();
		std::vector<double> IntBinMid; IntBinMid.resize(IntensityBins.size()); for(unsigned int i = 0; i < IntensityBins.size(); i++) IntBinMid.at(i) = (IntensityBins.at(i).first + IntensityBins.at(i).second)/2.;
		std::vector<double> Int; Int.resize(IntensityBins.size());
		std::vector<double> IntSum; IntSum.resize(IntensityBins.size());
		for(unsigned int i = 0; i < AstrophysicalSources.size(); i++)
		{							// calculate E^2*I / delta E  for every energy bin individually
			for(unsigned int j = 0; j < IntensityBins.size(); j++)	
			{
				Int.at(j) =pow(IntBinMid.at(j),2)/(IntensityBins.at(j).second - IntensityBins.at(j).first)* pow(10., Parameters.at(i)()) * AstrophysicalSources.at(i)->Intensity.at(j);
				IntSum.at(j) = ( i ==0 ? Int.at(j) : Int.at(j) + IntSum.at(j));
			}
			auto g = new TGraph(IntBinMid.size(), IntBinMid.data(), Int.data());
			g->SetName(("Int" + AstrophysicalSources.at(i)->Name).c_str()); g->Write();
		}
		int n = AstrophysicalSources.size();
		for(unsigned int i = 0; i < dmModels.size(); i++)
		{
			for(unsigned int j = 0; j < IntensityBins.size(); j++)
			{
				Int.at(j) = pow(IntBinMid.at(j),2)/(IntensityBins.at(j).second - IntensityBins.at(j).first) * pow(10, Parameters.at(i+n)()) * dmModels.at(i)->Intensity.at(j);
				IntSum.at(j) += Int.at(j);
			}			
			auto g = new TGraph(IntBinMid.size(), IntBinMid.data(), Int.data());
			g->SetName(("Int" + dmModels.at(i)->Name).c_str()); g->Write();
		}
		n += dmModels.size();
		for(unsigned int i = 0; i < AstrophysicalSourceClasses.size(); i++)
		{							// use the fitted second parameter
			for(unsigned int j = 0; j < IntensityBins.size(); j++)	
			{
				Int.at(j) = pow(IntBinMid.at(j),2)/(IntensityBins.at(j).second - IntensityBins.at(j).first)* pow(10., Parameters.at(n+2*i)()) * 
											AstrophysicalSourceClasses.at(i)->Intensity.at(j)->Eval(Parameters.at(2*i+n+1)());
				IntSum.at(j) += Int.at(j);
			}
			auto g = new TGraph(IntBinMid.size(), IntBinMid.data(), Int.data());
			g->SetName(("Int" + AstrophysicalSourceClasses.at(i)->Name).c_str()); g->Write();
		}
		auto gsum = new TGraph(IntBinMid.size(), IntBinMid.data(), IntSum.data());
		gsum->SetName("IntSum"); gsum->Write();
		/// now actual data
		std::vector<double> Intensity; Intensity.resize(IntensityBins.size()); 
		std::vector<double> IntensityErrors; IntensityErrors.resize(IntensityBins.size()); 
		std::vector<double> IntBinErrors; IntBinErrors.resize(IntensityBins.size());
		for(unsigned int i = 0; i < IntensityBins.size(); i++) 
		{
			Intensity.at(i) =  pow(IntBinMid.at(i),2)/(IntensityBins.at(i).second - IntensityBins.at(i).first) *IntensityMeasurement.at(i).first;
			IntensityErrors.at(i) = pow(IntBinMid.at(i),2)/(IntensityBins.at(i).second - IntensityBins.at(i).first) * IntensityMeasurement.at(i).second;
			IntBinErrors.at(i) = IntensityBins.at(i).second - IntensityBins.at(i).first;
		}
		auto gmeasure = new TGraphErrors(IntBinMid.size(), IntBinMid.data(), Intensity.data(), IntBinErrors.data(), IntensityErrors.data());
		gmeasure->SetName("Fermi"); gmeasure->Write();
	}
	
};


class IntensityAndAPSFit : public IntensityFit
{
protected:
	std::vector<Bounds> APSBins;
	std::shared_ptr<AngularPowerSpectrum<Measurement> > APSMeasurement;
	
	Bounds SBounds;
	
	double APSChiSquared()
	{
		double chiSquared = 0, val;
		//int nS = nPar; 	// position of S parameter
		double S_t = Parameters.at(Parameters.size()-1).getValue();	// S_t is last parameter
		//std::cout << "APSChiSquared:  S = " << S_t <<std::endl;
		for(unsigned int i = 0; i < APSMeasurement->Bin1Size(); i++)
		{
			for(unsigned int j = 0; j <= i; j++)
			{
				for(unsigned int k = 0; k < APSMeasurement->MultipoleNumber(); k++)
				{
					val = 0;
					for(unsigned int l =0; l < AstrophysicalSources.size(); l++)
					{
						val += pow(10, Parameters.at(l)())*AstrophysicalSources[l]->APS->at(i, j)->Eval(S_t);	// parameter is coefficient
					}
					int n = AstrophysicalSources.size();
					for(unsigned int l = 0; l < dmModels.size(); l++)
					{
						val += pow(10, 2.*Parameters.at(n+l)())* dmModels[l]->APS->at(i, j, k);		// this dependence is squared
					}
					n += dmModels.size();
					for(unsigned int l = 0; l < AstrophysicalSourceClasses.size(); l++)
					{
						auto asc = AstrophysicalSourceClasses[l];
						val += pow(10, Parameters.at(n+2*l)())*asc->APS->at(i, j)->Eval(S_t, Parameters.at(n+2*l+1)());	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
					}
					
					chiSquared += pow( (APSMeasurement->at(i, j, k).first - val)/APSMeasurement->at(i, j, k).second, 2);
				}
			}
		}
		int n = APSBins.size(); 	// number of elements in the triangular matrix should be given by (n^2 - n)/2
		return chiSquared/((pow(n,2) -n)/2*APSMeasurement->MultipoleNumber());
	}
	
	double loglike() override
	{
		return -0.5*(APSChiSquared() + IntensityChiSquared());
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
			
			Parameters.push_back(MultiNestParameter(SBounds, "S_t"));
			options.nlive = 500;
		} 
	
	void printResults()
	{
		std::cout << "IntensityAndAPS Fit: " << std::endl;
		for(unsigned int i = 0; i < Parameters.size(); i++) Parameters.at(i).print();
		std::cout << "chi^2 = " << IntensityChiSquared()+APSChiSquared() << std::endl;
	}
	
	void plotResults(TFile* f)
	{
		IntensityFit::plotResults(f);
		f->cd();
	}
	
};



#endif
