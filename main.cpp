#include "Benchmarking.h"
#include "LoadFromFiles.h"
#include "TApplication.h"
#include "TFile.h"
#include "TF1.h"
#include "CanvasWrapper.h"
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
	
	// Initialize CosmologyModel
	auto CM = std::make_shared<LambdaCDM>(Bounds(1e-3, 100), 100);
	
	// Load Linear Matter Power Spectrum
	std::vector<std::fstream> PkFiles;
	for(unsigned int i = 1; i <= 61; i++) 	// check file 62 and 63
		PkFiles.push_back(std::fstream(data + "LinearMatterPowerSpectrum/" + "honey_z" + std::to_string(i) + "_pk.dat", std::fstream::in));
	auto Plin = LoadLinearMatterPowerSpectrum(PkFiles);
	//Plin->print();
	
	// Load N of log(x) for DM
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
	auto DGRBAPS = LoadAnistropy(data + "Anisotropy2FGL/", 13,  20, Multipoles);
	
	// Prepare the Halo Model
	auto HM = std::make_shared<HaloModel>(CM, Plin, Bounds(1e-6_M_solar, 1e18_M_solar), Bounds(1e-3, 50.), Bounds(1e-2/1._Mpc, 1e3/1._Mpc));
	//HM->Init(30, 30, 50/*, data + "HM.dat"*/);
	HM->Load(data + "HM.dat");
	
	if(plot) plotHaloModel(HM, "plots/");
	//return 0;
	
	// Detector
	auto DT = std::make_shared<FermiLAT>();
	
	// Prepare Benchmark class
	Benchmark B(CM, HM, DT, true, plot, "plots/");
	B.LuminosityBounds_global = Bounds(1e30_ergpers, 1e52_ergpers); 	
	B.zBounds_global = Bounds(1e-3, 10);  B.zGridLen = 50;
	B.kBounds_global = Bounds(1e-2/1._Mpc, 9e2/1._Mpc); B.kGridLen = 200;
	B.SBounds_global = Bounds(1e-20_photonspercm2s, 1e-4_photonspercm2s);  B.SGridLen = 30;
	B.EBounds_global = Bounds(0.1_GeV, 900._GeV);  B.EGridLen = 50;
	B.GammaGridLen = 20;
	
	B.IntensityBins = IntensityBins;
	
	B.APSBins = APSBins;

	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::vector<std::shared_ptr<DarkMatter> > dmModels;
	
	AstrophysicalSources.push_back(std::make_shared<MAGN>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<FSRQ>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<LISP>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<HSP>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<NGSFG>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<SBSFG>(CM, tau));
	AstrophysicalSources.push_back(std::make_shared<SFGAGN>(CM, tau));
	B.calculateIntensityAndAPSForAstrophysicalSources(AstrophysicalSources); 
	
	//dmModels.push_back(std::make_shared<DecayingDM>(CM, HM, tau, dNdLogx, 10, 1.4e17));
	//dmModels.push_back(std::make_shared<AnnihilatingDM>(CM, HM, tau, dNdLogx, 10, 3e-26));
	
	//B.calculateIntensityAndAutocorrelationForDM(dmModels, Multipoles);
	
	return 0;
	
	std::vector<std::shared_ptr<AstrophysicalSourceClass> > ASC;	
	
	//auto IF = std::make_shared<IntensityFit>(IntensityBins, DGRBIntensity, Sources, ASC, "fit/ifit.dat") ;
	//IF->Run();
	//IF->printResults();
	
	auto IAF = std::make_shared<IntensityAndAPSFit>(IntensityBins, DGRBIntensity, APSBins, DGRBAPS, AstrophysicalSources, dmModels,  ASC, Bounds(1e-10, 4e-10), "fit/iafit.dat");
	IAF->Run();
	IAF->printResults();
	return 0;
}

void plotHaloModel(std::shared_ptr<HaloModel> HM, std::string path)
{
	auto rootFile = new TFile((path + "HaloModelPlots.root").c_str(), "RECREATE", "HaloModelPlots");
	
	/// Linear Halo Bias
	Canvas LinearHaloBias("LHB", "Linear Halo Bias Plot", 1000, 1000);
	LinearHaloBias().SetLogx(1); LinearHaloBias().SetLogy(1);
	double* z = new double[4]; z[0]= 0; z[1] = 1; z[2] = 2; z[3] = 3;
	
	for(unsigned int i = 0; i < 4; i++)
	{
		TF1 LHB(("LHB2z" + std::to_string(z[i])).c_str(), [HM, z, i] (double* args, double* params) { return pow(HM->LinearHaloBias(args[0]*h, z[i]),2.); },
						1e10, 1e14, 0);
		LHB.SetNpx(1e5);
		auto g = new TGraph(&LHB);  //g->SetLineColor(i+1);
		LinearHaloBias.AddGraph(g, ("z = " + std::to_string(int(z[i]))).c_str(), "L", true);
	}
	LinearHaloBias.Draw("A");	
	LinearHaloBias.SetxLimits(1e10, 1e14);
	LinearHaloBias.SetyLimits(1e-1, 1e2);
	LinearHaloBias().BuildLegend();
	LinearHaloBias().SaveAs((path  + "LHB.jpg").c_str());
	
	/// Halo Mass Function
	Canvas HMFCanvas("HMF", "Halo Mass Function Plot", 1000, 1000);
	HMFCanvas().SetLogx(1); HMFCanvas().SetLogy(1);
	
	unsigned int n = 30;  double M_min = 1e4;  double M_max = 1e12;
	double* M = new double[n];
	for(unsigned int i =0; i < n; i++) M[i] = exp( log(M_min) + i*(log(M_max) - log(M_min))/(n-1.));
	double* HMF = new double[n];
	z[0] = 0; z[1] =  10; z[2] = 20; z[3] =  30;
	for(unsigned int i = 0; i < 4; i++)
	{
		for(unsigned int j = 0; j < n; j++) HMF[j] = HM->HaloMassFunction(M[j]*h, z[i]);
		auto HMFGraph = new TGraph(n, M, HMF);
		//HMFGraph->SetTitle(("z = " + std::to_string(z[i])).c_str());  
		//HMFGraph->SetLineColor(i+1);
		//HMFGraph->Draw((std::string(i ==0 ? "A" : "same") + "L").c_str());
		//auto g = new TGraph(&HMF);  g->SetLineColor(i+1);
		HMFCanvas.AddGraph(HMFGraph, ("z = " + std::to_string(int(z[i]))).c_str(), "L", true);
	}
	
	HMFCanvas.Draw("A");	
	HMFCanvas.SetxLimits(M_min, M_max);
	HMFCanvas.SetyLimits(1e-7, 1e5);
	HMFCanvas.SetAxesTitle("M", "dn/dM");
	HMFCanvas().BuildLegend();
	HMFCanvas().SaveAs((path + "HMF.jpg").c_str());
	
	delete []M; delete []HMF;
	delete []z;
	
	/// Integrated Halo Mass Function  
	/*Canvas IHMFCanvas("nCanvas", "Integrated Halo Mass Function", 1000, 1000);
	IHMFCanvas().SetLogy(1);
	n = 40;  z = new double[n]; for(unsigned int i = 0; i < n; i++) z[i] = i*(40.)/(n-1.);
	int* Mthres = new int[11]; for(unsigned int i = 0; i < 11; i++) Mthres[i] = 4 + i;
	double* IHMF = new double[n];
	
	TF1 _HMF("HMF", [HM] (double* args, double* params) { return exp(args[0])*HM->HaloMassFunction(h*exp(args[0]), params[0]); }, log(1e4), log(1e18), 1);
	
	for(unsigned int i = 0; i < 11; i++)
	{
		
		for(unsigned int j =0; j < n; j++) 
		{
			_HMF.SetParameters(z[j], 0.);
			IHMF[j] = _HMF.Integral(log(pow(10., Mthres[i])), log(1e18), 1e-4);
		}
		auto IHMFGraph = new TGraph(n, z, IHMF);  IHMFGraph->SetLineColor(i+1);
		IHMFCanvas.AddGraph(IHMFGraph, "$M_{thres}$ = 1e" + std::to_string(Mthres[i]), "L");
	}
	IHMFCanvas.Draw("A");
	IHMFCanvas.SetyLimits(1e-10, 1e10);
	IHMFCanvas.SetxLimits(-1, 40);
	IHMFCanvas().BuildLegend();
	IHMFCanvas().SaveAs("IHMF.jpg");
	
	delete []Mthres; delete []z; */
	
	rootFile->Close();
	delete rootFile;	
}

void plotIntensities(std::vector<std::shared_ptr<DGRBSource> > sources, std::vector<Bounds> IntensityBins, std::string path)
{
	std::vector<double> EBinsMid; EBinsMid.resize(IntensityBins.size());
	for(unsigned int i =0; i < IntensityBins.size(); i++) EBinsMid[i] = (IntensityBins[i].first + IntensityBins[i].second)/2.;
	std::vector<double> Intensities; Intensities.resize(IntensityBins.size());
	Canvas IntensityCanvas("Int", "Intensity");
	IntensityCanvas().SetLogx(1); IntensityCanvas().SetLogy(1); 
	for(unsigned int i = 0; i < sources.size(); i++)
	{
		for(unsigned int j = 0; j < IntensityBins.size(); j++) Intensities[j] = sources[i]->Intensity[j]*pow(EBinsMid[j], 2.)/(IntensityBins[j].second - IntensityBins[j].first); 
		TGraph* g = new TGraph(IntensityBins.size(), EBinsMid.data(), Intensities.data());
		IntensityCanvas.AddGraph(g, sources[i]->Name, "L", true);
	}
	IntensityCanvas.Draw("A");
	IntensityCanvas.SetAxesTitle("E / GeV", "I*$E^2$/$\\Delta$E");
	IntensityCanvas().BuildLegend();
	IntensityCanvas().SaveAs((path + "ASIntensity.jpg").c_str());
}


