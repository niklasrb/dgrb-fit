#include "Benchmarking.h"
#include "LoadFromFiles.h"
#include "TFile.h"
#include "dgrbFit.h"




// g++ -Wall -c "%f" -std=c++11 -I$ROOTSYS/include -g
// g++ -Wall -v -o "%e" "%f" $(sh $ROOTSYS/bin/root-config --cflags --glibs) -lgsl -lgslcblas -L/usr/local/lib/libnest/ -ggdb -lgfortran -lnest3 -llapack 

void plotHaloModel(std::shared_ptr<HaloModel> HM, std::string path);
void plotIntensities(std::vector<std::shared_ptr<DGRBSource> > sources, std::vector<Bounds> IntensityBins, std::string path);


int main(int argc, char** argv)
{
	bool plot = true;
	
	std::string data("data/");
	
	// Load model values for the extra galactic background absorbtion
	std::fstream dominguez(data + std::string("dominguez11.txt"),  std::fstream::in);
	auto tau = LoadEBLAbsorbtionCoefficient(dominguez);
	if(plot) tau->plot("plots/EBLAbsCoeff.root");
	
	// Initialize CosmologyModel
	auto CM = std::make_shared<LambdaCDM>(Bounds(1e-3, 100), 100);
	
	// Load Linear Matter Power Spectrum
	std::vector<std::fstream> PkFiles;
	for(unsigned int i = 1; i <= 61; i++) 	// check file 62 and 63
		PkFiles.push_back(std::fstream(data + "LinearMatterPowerSpectrum/" + "honey_z" + std::to_string(i) + "_pk.dat", std::fstream::in));
	auto Plin = LoadLinearMatterPowerSpectrum(PkFiles);
	if(plot) Plin->plot("plots/LMPS.root");
	
	// Load dN/dlogX for DM
	std::fstream dmEnergySpectrum(data + "AtProduction_gammas.dat", std::fstream::in);
	auto dNdLogx = LoaddNdLogx(dmEnergySpectrum);
	//dNdLogx->print();
	
	// Load DGRB intensity data
	std::vector<Measurement> DGRBIntensity;
	std::vector<Bounds> IntensityBins;
	std::fstream dgrbIntensityFile(data + "DGRBIntensity.txt", std::fstream::in);
	LoadDGRB(dgrbIntensityFile, IntensityBins, DGRBIntensity);
	
	// Load DGRB anisotropy
	std::vector<double> Multipoles;
	std::vector<Bounds> APSBins = {Bounds(0.5_GeV, 0.72_GeV ), Bounds(0.72_GeV, 1.04_GeV), Bounds(1.04_GeV, 1.38_GeV), Bounds(1.38_GeV, 1.99_GeV), 
									Bounds(1.99_GeV, 3.15_GeV), Bounds(3.15_GeV, 5.0_GeV), Bounds(5.0_GeV, 7.23_GeV), Bounds(7.23_GeV, 10.24_GeV), 
									Bounds(10.24_GeV, 21.83_GeV), Bounds(21.83_GeV, 50.0_GeV), Bounds(50.0_GeV, 95.27_GeV), Bounds(95.27_GeV, 199.05_GeV),
									Bounds(199.05_GeV, 500.0_GeV)  };
	auto _2FGLAPS = LoadAnistropy(data + "Anisotropy2FGL/", 13,  20, Multipoles);	// these are the whole Cl data
	auto _2FGLCp = LoadAnistropy(data + "Anisotropy2FGL/Cp.txt", 13);				// this is Cp data
	
	// Prepare the Halo Model
	auto HM = std::make_shared<HaloModel>(CM, Plin, Bounds(1e-6_M_solar, 1e18_M_solar), Bounds(9e-4, 50.), Bounds(9e-4/1._Mpc, 2e4/1._Mpc));
	HM->Init(50, 50, 100, data + "HM.dat");			// recalculate from scratch and save in file
	//HM->Load(data + "HM.dat");							// load from file
	
	if(plot) HM->Plot("plots/HM.root");
	
	// Detector
	auto DT = std::make_shared<FermiLAT>();
	
	// Prepare Benchmark class
	Benchmark B(CM, HM, DT, true, plot, "plots/");
	B.LuminosityBounds_global = Bounds(1e20_ergpers, 1e52_ergpers); 	
	B.zBounds_global = Bounds(1e-3, 10);  B.zGridLen = 50;
	B.kBounds_global = Bounds(1e-3/1._Mpc, 1e4/1._Mpc); B.kGridLen = 200;
	B.SBounds_global = Bounds(1e-24_photonspercm2s, 1e-4_photonspercm2s);  B.SGridLen = 50;
	B.EBounds_global = Bounds(0.1_GeV, 900._GeV);  B.EGridLen = 50;
	B.GammaGridLen = 30;
	
	B.IntensityBins = IntensityBins;
	
	B.APSBins = APSBins;

	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::vector<std::shared_ptr<DarkMatter> > dmModels;
	std::vector<std::shared_ptr<AstrophysicalSourceClass> > ASC;	
	
	/// add sources
	AstrophysicalSources.push_back(std::make_shared<MAGN>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<FSRQ>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<LISP>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<HSP>(CM, tau));
	auto SB = std::make_shared<SBSFG>(CM, tau); AstrophysicalSources.push_back(SB);
	auto NG = std::make_shared<NGSFG>(CM, tau); AstrophysicalSources.push_back(NG);
	auto AGN = std::make_shared<SFGAGN>(CM, tau); AstrophysicalSources.push_back(AGN);
	/// calculate intensity and anisotropy for them
	B.calculateIntensityAndAPSForAstrophysicalSources(AstrophysicalSources);
	
	return 0;
	
	AstrophysicalSources.erase(AstrophysicalSources.begin() + 1 , AstrophysicalSources.end());	// delete SFG elements
	AstrophysicalSources.push_back(std::make_shared<SFG>(SB, NG, AGN));						// and add combination
	 
	/// test dependency on E_cut for FSRQ
	std::vector<double> FSRQE_cut = {2.5, 12.5, 22.5, 32.5, 42.5, 52.5, 62.5};
	std::vector<std::shared_ptr<AstrophysicalSource> > FSRQs;
	for(unsigned int i = 0; i < FSRQE_cut.size(); i++) { FSRQs.push_back(std::make_shared<FSRQ>(CM, tau, FSRQE_cut.at(i))); FSRQs.at(i)->Name += std::to_string(FSRQE_cut.at(i)); }
	B.calculateIntensityAndAPSForAstrophysicalSources(FSRQs);
	
	ASC.push_back(std::make_shared<AstrophysicalSourceClass>(IntensityBins, APSBins, FSRQs, FSRQE_cut, "FSRQ"));
	//AstrophysicalSources.erase(AstrophysicalSources.begin() + 1);
	
	/// test dependency on E_cut for HSP
	std::vector<double> HSPE_cut = {100., 300., 500., 700., 900., 1100., 1500., 2000., 2600.};
	std::vector<std::shared_ptr<AstrophysicalSource> > HSPs;
	for(unsigned int i = 0; i < HSPE_cut.size(); i++) { HSPs.push_back(std::make_shared<HSP>(CM, tau, HSPE_cut.at(i))); HSPs.at(i)->Name += std::to_string(HSPE_cut.at(i)); }
	B.calculateIntensityAndAPSForAstrophysicalSources(HSPs);
	
	ASC.push_back(std::make_shared<AstrophysicalSourceClass>(IntensityBins, APSBins, HSPs, HSPE_cut, "HSP"));
	//AstrophysicalSources.erase(AstrophysicalSources.begin() + 1);
	
	/// test dependency on E_cut for LISP
	std::vector<double> LISPE_cut = {20, 30, 40, 50, 60, 70, 80};
	std::vector<std::shared_ptr<AstrophysicalSource> > LISPs;
	for(unsigned int i = 0; i < LISPE_cut.size(); i++) { LISPs.push_back(std::make_shared<LISP>(CM, tau, LISPE_cut.at(i))); LISPs.at(i)->Name += std::to_string(LISPE_cut.at(i)); }
	B.calculateIntensityAndAPSForAstrophysicalSources(LISPs);
	
	ASC.push_back(std::make_shared<AstrophysicalSourceClass>(IntensityBins, APSBins, LISPs, LISPE_cut, "LISP"));
	//AstrophysicalSources.erase(AstrophysicalSources.begin() + 1);
	
	
	/// add dark matter models
	//dmModels.push_back(std::make_shared<DecayingDM>(CM, HM, tau, dNdLogx, 10, 1.4e17)); dmModels.at(0)->Name = "M=10GeV";
	//dmModels.push_back(std::make_shared<DecayingDM>(CM, HM, tau, dNdLogx, 100, 1.4e17)); dmModels.at(1)->Name = "M=100GeV";
	//dmModels.push_back(std::make_shared<DecayingDM>(CM, HM, tau, dNdLogx, 1000, 1.4e17)); dmModels.at(2)->Name = "M=1TeV";
	
	//dmModels.push_back(std::make_shared<AnnihilatingDM>(CM, HM, tau, dNdLogx, 10, 3e-26));
	
	//B.calculateIntensityAndAutocorrelationForDM(dmModels, Multipoles);
	
	//return 0;	
	/*
	auto IF = std::make_shared<IntensityFit>(IntensityBins, DGRBIntensity, AstrophysicalSources, dmModels, ASC, "fits/intfit/intfitase") ;
	IF->Run();
	IF->printResults();
	auto f = new TFile("plots/intfitase.root", "RECREATE"); IF->plotResults(f); delete f;
	*/
	/*
	auto IAF = std::make_shared<IntensityAndAPSFit>(IntensityBins, DGRBIntensity, APSBins, _2FGLCp, std::vector<double>(), AstrophysicalSources, dmModels,  ASC, Bounds(3e-10, 6e-10), "fits/intapsfit/as");
	IAF->Run();
	IAF->printResults();
	auto f = new TFile("plots/apsfit.root", "RECREATE"); IAF->plotResults(f); delete f;
	*/
	/*
	auto IAPSF = std::make_shared<IntensityAndAPSFit>(IntensityBins, DGRBIntensity, APSBins, _2FGLAPS, Multipoles, AstrophysicalSources, dmModels,  ASC, Bounds(3e-10, 6e-10), "fits/intapsfit/all");
	IAPSF->Run();
	IAPSF->printResults();
	auto f2 = new TFile("plots/apsfitall.root", "RECREATE"); IAPSF->plotResults(f2); delete f2;
	*/
	return 0;
}




