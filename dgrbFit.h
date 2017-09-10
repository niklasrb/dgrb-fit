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
	
	double IntensityChiSquared();
	double loglike() override;
	void dumper(MultiNestFitData& data) override;	
	
public:
	IntensityFit(const std::vector<Bounds>& IntensityBins, const std::vector<Measurement>& IntensityMeasurement, const std::vector<std::shared_ptr<AstrophysicalSource> >& AstrophysicalSources, const std::vector<std::shared_ptr<DarkMatter> >& dmModels,
					const std::vector<std::shared_ptr<AstrophysicalSourceClass> >& AstrophysicalSourceClasses, std::string output) ;
	
	void printResults();
	void plotResults(TFile* f);

};

class IntensityAndAPSFit : public IntensityFit
{
protected:
	std::vector<Bounds> APSBins;
	std::shared_ptr<AngularPowerSpectrum<Measurement> > APSMeasurement;
	
	Bounds SBounds;
	
	// these functions expect pointers, so that they can be used if there are two APS measurements in the class, i.e. 2FGL and 3FGL catalogue
	double APSChiSquared(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t);		// handles multipole depenedance
	double CpChiSquared(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t);		// only use this if there is no multipole dependance // ignores DM!
	double loglike() override;
	void dumper(MultiNestFitData& data) override;
	void plotCp(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t, TFile* f);
	void plotAPS(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t, TFile* f);
	
public:
	IntensityAndAPSFit(const std::vector<Bounds>& IntensityBins, const std::vector<Measurement>& IntensityMeasurement, const std::vector<Bounds>& APSBins, const std::shared_ptr<AngularPowerSpectrum<Measurement> >& APSMeasurement,
	const std::vector<std::shared_ptr<AstrophysicalSource> >& AstrophysicalSources, const std::vector<std::shared_ptr<DarkMatter> >& dmModels, const std::vector<std::shared_ptr<AstrophysicalSourceClass> >& AstrophysicalSourceClasses, Bounds SBounds, std::string output);

	void printResults();
	void plotResults(TFile* f);
};

/// Intensity Fit

double IntensityFit::IntensityChiSquared()
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
	
double IntensityFit::loglike()
{
	return -0.5*IntensityChiSquared();
}
	
void IntensityFit::dumper(MultiNestFitData& data)
{
	
	std::cout << "Last parameters: "; for(int i = 0; i < data.nPar; i++) std::cout << i << ": "<< data.lastParameters[i] << (i < data.nPar-1 ? '\t' : '\n');
	std::cout << "loglike: " << data.logLike << std::endl;
}


IntensityFit::IntensityFit(const std::vector<Bounds>& IntensityBins, const std::vector<Measurement>& IntensityMeasurement, const std::vector<std::shared_ptr<AstrophysicalSource> >& AstrophysicalSources, const std::vector<std::shared_ptr<DarkMatter> >& dmModels,
				const std::vector<std::shared_ptr<AstrophysicalSourceClass> >& AstrophysicalSourceClasses, std::string output = "") 
: MultiNestFit(output) , IntensityBins(IntensityBins), IntensityMeasurement(IntensityMeasurement), AstrophysicalSources(AstrophysicalSources), dmModels(dmModels), AstrophysicalSourceClasses(AstrophysicalSourceClasses)
{
	assert(IntensityBins.size() == IntensityMeasurement.size() && IntensityBins.size() > 0);
	for(unsigned int i = 0; i < AstrophysicalSources.size(); i++) Parameters.push_back(MultiNestParameter(Bounds(-2, 3), AstrophysicalSources[i]->Name + " LumFunc log factor"));
	for(unsigned int i = 0; i < dmModels.size(); i++) Parameters.push_back(MultiNestParameter(Bounds(-5, 5), dmModels[i]->Name + " Windowfunction log factor"));
	for(unsigned int i = 0; i < AstrophysicalSourceClasses.size(); i++)
	{
		Parameters.push_back(MultiNestParameter(Bounds(-2, 3), AstrophysicalSourceClasses[i]->Name + " LumFunc log factor"));
		Parameters.push_back(MultiNestParameter(AstrophysicalSourceClasses[i]->parameterBounds, AstrophysicalSourceClasses[i]->Name + " second param"));
	}
}
	
void IntensityFit::printResults()
{
	std::cout << "Intensity Fit: " << std::endl;
	for(unsigned int i = 0; i < Parameters.size(); i++) Parameters.at(i).print();
	std::cout << "chi^2 = " << IntensityChiSquared() << std::endl;
}

void IntensityFit::plotResults(TFile* f)
{
	f->cd();
	std::vector<double> IntBinMid; IntBinMid.resize(IntensityBins.size()); 
	for(unsigned int i = 0; i < IntensityBins.size(); i++) IntBinMid.at(i) = (IntensityBins.at(i).first + IntensityBins.at(i).second)/2.;
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
		IntBinErrors.at(i) = (IntensityBins.at(i).second - IntensityBins.at(i).first)/2.;
		IntensityErrors.at(i) = pow(IntBinMid.at(i),2)/(IntBinErrors.at(i)*2) * IntensityMeasurement.at(i).second;
	}
	auto gmeasure = new TGraphErrors(IntBinMid.size(), IntBinMid.data(), Intensity.data(), IntBinErrors.data(), IntensityErrors.data());
	gmeasure->SetName("FermiInt"); gmeasure->Write();
}
	
/// IntensityAndAPS Fit

/// This function calculates (sum - measurement)^2 / (measurement_error)  for all energy bins i,j with j <=i and for every multipole
/// 	and adds it together to get chi^2, then dividing that by the number of energy bin combinations*multipole
double IntensityAndAPSFit::APSChiSquared(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t)
{
	double chiSquared = 0, sum;
	for(unsigned int i = 0; i < APSMeasurement->Bin1Size(); i++)
	{
		for(unsigned int j = 0; j <= i; j++)
		{
			for(unsigned int k = 0; k < APSMeasurement->MultipoleNumber(); k++)
			{
				sum = 0;
				for(unsigned int l =0; l < AstrophysicalSources.size(); l++)
				{
					sum += pow(10, Parameters.at(l)())*AstrophysicalSources.at(l)->APS->at(i, j)->Eval(S_t);	// parameter is coefficient
				}
				int n = AstrophysicalSources.size();
				for(unsigned int l = 0; l < dmModels.size(); l++)
				{
					sum += pow(10, 2.*Parameters.at(n+l)())* dmModels.at(l)->APS->at(i, j, k);		// this dependence is squared
				}
				n += dmModels.size();
				for(unsigned int l = 0; l < AstrophysicalSourceClasses.size(); l++)
				{
					auto asc = AstrophysicalSourceClasses.at(l);
					sum += pow(10, Parameters.at(n+2*l)())*asc->APS->at(i, j)->Eval(S_t, Parameters.at(n+2*l+1)());	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
				}
				
				chiSquared += pow( (APS->at(i, j, k).first - sum)/APS->at(i, j, k).second, 2);
			}
		}
	}
	int n = APSBins.size(); 	// number of elements in the triangular matrix should be given by (n^2 - n)/2
	return chiSquared/((pow(n,2) -n)/2*APSMeasurement->MultipoleNumber());
}

/// This function calculates (sum - measurement)^2 / (measurement_error)  for all energy bins i,j with j <=i 
///		and adds it together to get chi^2, then dividing that by the number of energy bin combinations
/// ignores dark matter components
double IntensityAndAPSFit::CpChiSquared(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t)
{
	double chiSquared = 0, sum;
	for(unsigned int i = 0; i < APSMeasurement->Bin1Size(); i++)
	{
		for(unsigned int j = 0; j <= i; j++)
		{
			sum = 0;
			for(unsigned int l =0; l < AstrophysicalSources.size(); l++)
			{
				sum += pow(10, Parameters.at(l)())*AstrophysicalSources.at(l)->APS->at(i, j)->Eval(S_t);	// parameter is coefficient
			}
			int n = AstrophysicalSources.size() + dmModels.size();
			for(unsigned int l = 0; l < AstrophysicalSourceClasses.size(); l++)
			{
				auto asc = AstrophysicalSourceClasses.at(l);
				sum += pow(10, Parameters.at(n+2*l)())*asc->APS->at(i, j)->Eval(S_t, Parameters.at(n+2*l+1)());	// first parameter is coefficient, second parameter is class parameter, i.e. E_cut
			}
			chiSquared += pow( (APS->at(i, j, 0).first - sum)/APS->at(i, j, 0).second, 2);
		}
	}
	int n = APSBins.size(); 	// number of elements in the triangular matrix should be given by (n^2 - n)/2
	return chiSquared/((pow(n,2) -n)/2.);
}

double IntensityAndAPSFit::loglike()
{
	const double aps = ( APSMeasurement->MultipoleNumber() == 1 ? CpChiSquared(APSMeasurement, Parameters.at(Parameters.size()-1)()) 
																: APSChiSquared(APSMeasurement, Parameters.at(Parameters.size()-1)()));	// differentiate wether there is multipole dependance
	return -0.5*(aps + IntensityChiSquared());
}

void IntensityAndAPSFit::dumper(MultiNestFitData& data)
{
	
	std::cout << "Last parameters: "; for(int i = 0; i < data.nPar; i++) std::cout << i << ": "<< data.lastParameters[i] << (i < data.nPar-1 ? '\t' : '\n');
	std::cout << "loglike: " << data.logLike << std::endl;
}

IntensityAndAPSFit::IntensityAndAPSFit(const std::vector<Bounds>& IntensityBins, const std::vector<Measurement>& IntensityMeasurement, const std::vector<Bounds>& APSBins, const std::shared_ptr<AngularPowerSpectrum<Measurement> >& APSMeasurement,
	const std::vector<std::shared_ptr<AstrophysicalSource> >& AstrophysicalSources, const std::vector<std::shared_ptr<DarkMatter> >& dmModels, const std::vector<std::shared_ptr<AstrophysicalSourceClass> >& AstrophysicalSourceClasses, Bounds SBounds, std::string output = "")
	: IntensityFit(IntensityBins, IntensityMeasurement, AstrophysicalSources, dmModels, AstrophysicalSourceClasses, output) , APSBins(APSBins), APSMeasurement(APSMeasurement), SBounds(SBounds)
{
	assert(APSBins.size() == APSMeasurement->Bin1Size()  && APSBins.size() == APSMeasurement->Bin2Size() && APSBins.size() >=2);
	
	Parameters.push_back(MultiNestParameter(SBounds, "S_t"));
	options.nlive = 1500;
} 

void IntensityAndAPSFit::printResults()
{
	std::cout << "IntensityAndAPS Fit: " << std::endl;
	for(unsigned int i = 0; i < Parameters.size(); i++) Parameters.at(i).print();
	std::cout << "chi^2 = " << -2.*loglike() << std::endl;
}

void IntensityAndAPSFit::plotResults(TFile* f)
{
	IntensityFit::plotResults(f);
	if(APSMeasurement->MultipoleNumber() == 1)	
		plotCp(APSMeasurement, Parameters.at(Parameters.size()-1)(), f);
	else
		plotAPS(APSMeasurement, Parameters.at(Parameters.size()-1)(), f);
}

void IntensityAndAPSFit::plotAPS(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t, TFile* f)
{
	/// plot shit for each energy bin
}

/// plots the Cp*E^4/DelaE^2 depending on the energy bin for all the sources and the observation data
void IntensityAndAPSFit::plotCp(std::shared_ptr<AngularPowerSpectrum<Measurement>> APS, const double& S_t, TFile* f)
{
	std::vector<double> APSBinMid, APSBinSize, APSBinSizeHalfed; APSBinMid.resize(APSBins.size()); APSBinSize.resize(APSBins.size()); APSBinSizeHalfed.resize(APSBins.size());
	std::vector<double> Cp, Cp_err, CpSum; Cp.resize(APSBins.size()); Cp_err.resize(APSBins.size()); CpSum.resize(APSBins.size());
	for(unsigned int i = 0; i < APSBins.size(); i++)
	{
		APSBinMid.at(i) = (APSBins.at(i).first + APSBins.at(i).second)/2.;
		APSBinSize.at(i) = APSBins.at(i).second - APSBins.at(i).first;
		APSBinSizeHalfed.at(i) = APSBinSize.at(i)/2.;
		Cp.at(i) = pow(APSBinMid.at(i), 4)/pow(APSBinSize.at(i), 2) * APS->at(i, i, 0).first;  
		Cp_err.at(i) = pow(APSBinMid.at(i), 4)/pow(APSBinSize.at(i), 2) * APS->at(i, i, 0).second;
		CpSum.at(i) = 0;
	}
	f->cd();
	auto gMeas = new TGraphErrors(APSBins.size(), APSBinMid.data(), Cp.data(), APSBinSizeHalfed.data(), Cp_err.data());
	gMeas->SetName("FermiCp"); gMeas->Write();
	// take Cp for astrophysical sources
 	for(unsigned int i = 0; i < AstrophysicalSources.size(); i++)
	{
		for(unsigned int j = 0; j < APSBins.size(); j++)	
			{	Cp.at(j) = pow(APSBinMid.at(j), 4)/pow(APSBinSize.at(j), 2) * pow(10, Parameters.at(i)())*AstrophysicalSources.at(i)->APS->at(j, j)->Eval(S_t); CpSum.at(j) += Cp.at(j); }
		auto g = new TGraph(APSBins.size(), APSBinMid.data(), Cp.data());
		g->SetName(("Cp" + AstrophysicalSources.at(i)->Name).c_str()); g->Write();
	}
	// plot astrophysical source classes
	int n = AstrophysicalSources.size() + dmModels.size();
	for(unsigned int i = 0; i < AstrophysicalSourceClasses.size(); i++)
	{
		auto asc = AstrophysicalSourceClasses.at(i);
		for(unsigned int j = 0; j < APSBins.size(); j++)	
			{	Cp.at(j) = pow(APSBinMid.at(j), 4)/pow(APSBinSize.at(j), 2) * pow(10, Parameters.at(n+2*i)())*asc->APS->at(j, j)->Eval(S_t, Parameters.at(n+2*i+1)()); CpSum.at(j) += Cp.at(j); }
		auto g = new TGraph(APSBins.size(), APSBinMid.data(), Cp.data());
		g->SetName(("Cp" + asc->Name).c_str()); g->Write();
	}
	// plot sum
	auto g = new TGraph(APSBins.size(), APSBinMid.data(), CpSum.data());
	g->SetName("CpSum"); g->Write();
}




#endif
