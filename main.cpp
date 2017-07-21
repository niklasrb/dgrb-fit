#include "Benchmarking.h"
#include "LoadFromFiles.h"
#include "TApplication.h"
#include "TFile.h"
//#include "dgrbFit.h"

int main(int argc, char** argv)
{
	bool plot =true;
	//TFile* plotFile;
	//if(plot) plotFile = new TFile("plots.root", "RECREATE");
	TApplication* dummyApp;
	if(plot) dummyApp = new TApplication("dummy", &argc, argv);
	
	std::string data("/media/data/Documents/University/BachelorThesis/trafunc/");
	
	// Load model values for the extra galactic background absorbtion
	std::fstream dominguez(data + std::string("dominguez11.txt"),  std::fstream::in);
	auto tau = LoadEBLAbsorbtionCoefficient(dominguez);
	//tau->print();
	
	auto CM = std::make_shared<LambdaCDM>();
	
	// Load Linear Matter Power Spectrum
	std::vector<std::fstream> PkFiles;
	for(unsigned int i = 1; i <= 61; i++) 	// check file 62 and 63
		PkFiles.push_back(std::fstream(data + std::string("honey_z") + std::to_string(i) + std::string("_pk.dat"), std::fstream::in));
	auto Plin = LoadLinearMatterPowerSpectrum(PkFiles);
	//Plin->print();
	
	// Prepare the Halo Model
	auto HM = std::make_shared<HaloModel>(CM, Plin, Bounds(1e-6 * M_solar, 1e18 * M_solar), Bounds(1e-3, 6.));
	HM->Init(60, 13);
	
	// Define energy bins
	std::vector<Bounds> EBins = {Bounds(0.1_GeV ,.14_GeV), Bounds(0.14_GeV ,0.2_GeV), Bounds(0.2_GeV ,0.28_GeV), Bounds(0.28_GeV ,0.4_GeV),
								Bounds(0.4_GeV ,0.57_GeV), Bounds(0.57_GeV ,0.8_GeV),  Bounds(0.8_GeV ,1.1_GeV), Bounds(1.1_GeV ,1.6_GeV), 
								Bounds(1.6_GeV, 2.3_GeV), Bounds(2.3_GeV, 3.2_GeV), Bounds(3.2_GeV, 4.5_GeV), Bounds(4.5_GeV, 6.4_GeV),
								Bounds(6.4_GeV, 9.1_GeV), Bounds(9.1_GeV, 13._GeV), Bounds(13._GeV, 18._GeV), Bounds(18._GeV, 26._GeV),
								Bounds(26._GeV, 36._GeV), Bounds(36._GeV, 51._GeV), Bounds(51._GeV, 72._GeV), Bounds(72._GeV, 100._GeV),
								Bounds(100._GeV, 140._GeV), Bounds(140._GeV, 200._GeV), Bounds(200._GeV, 290._GeV), Bounds(290._GeV, 410._GeV),
								Bounds(410._GeV, 580._GeV), Bounds(580._GeV, 820._GeV)};
	
	// Prepare Benchmark class
	
	Benchmark B(CM, HM, EBins, Bounds(1e-3, 1e15), Bounds(1e-3, 100), true);
	//auto magn = std::make_shared<MAGN>(CM);
	//std::cout << magn->RescaledLuminosityFunction(1e5, 1, 2.37) << std::endl;
	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::vector<std::shared_ptr<DarkMatter> > dmModels;
	
	AstrophysicalSources.push_back(std::make_shared<MAGN>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<FSRQ>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<LISP>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<HSP>(CM, tau));
	//AstrophysicalSources.push_back(std::make_shared<SFG>(CM, tau));
	B.calculateIntensityAndAutocorrelationForAstrophysicalSources(AstrophysicalSources, 30, 8, false);
	
	//dmModels.push_back(std::make_shared<AnnihilatingDM>(CM, HM, tau, 1, 1));
	//dmModels.push_back(std::make_shared<DecayingDM>(CM, HM, tau, 1, 1));
	
	//B.calculateIntensityAndAutocorrelationForDM(dmModels, 20);
	
	for(unsigned int i = 0; i < AstrophysicalSources.size(); i++)
	{
		std::cout << AstrophysicalSources.at(i)->Name << ": " << std::endl;
		for(unsigned int j = 0; j < EBins.size(); j++) std::cout << AstrophysicalSources.at(i)->Intensity.at(j).second << '\t';
		std::cout << std::endl;
	}
	if(plot) 
	{
		//plotFile->Write();
		//delete plotFile;
		dummyApp->Run(true);
	}
	return 0;
}

