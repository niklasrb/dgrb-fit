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
void plotIntensities(std::vector<std::shared_ptr<AstrophysicalSource> >, std::string path);


int main(int argc, char** argv)
{
	bool plot = true;
	
	std::string data("data/");
	
	// Load model values for the extra galactic background absorbtion
	std::fstream dominguez(data + std::string("dominguez11.txt"),  std::fstream::in);
	auto tau = LoadEBLAbsorbtionCoefficient(dominguez);
	
	// Initialize CosmologyModel
	auto CM = std::make_shared<LambdaCDM>(Bounds(1e-3, 100), 100);
	//std::cout << "Critical Density: " << CM->CriticalDensity << std::endl;
	
	// Load Linear Matter Power Spectrum
	std::vector<std::fstream> PkFiles;
	for(unsigned int i = 1; i <= 61; i++) 	// check file 62 and 63
		PkFiles.push_back(std::fstream(data + "honey_z" + std::to_string(i) + "_pk.dat", std::fstream::in));
	auto Plin = LoadLinearMatterPowerSpectrum(PkFiles);
	//Plin->print();
	
	// Load N of log(x) for DM
	std::fstream dmEnergySpectrum(data + "AtProduction_gammas.dat", std::fstream::in);
	auto dNdLogx = LoaddNdLogx(dmEnergySpectrum);
	//dNdLogx->print();
	
	
	// Prepare the Halo Model
	auto HM = std::make_shared<HaloModel>(CM, Plin, Bounds(1e-6 /* solar masses*/, 1e18), Bounds(1e-3, 50.), Bounds(1e-2, 1e3));
	HM->Init(500, 500, 500, data + "HM.dat");
	
	//if(plot) plotHaloModel(HM, "");
	return 0;
	
	// Detector
	auto DT = std::make_shared<FermiLAT>();
	
	// Prepare Benchmark class
	Benchmark B(CM, HM, DT, true, plot);
	B.LuminosityBounds_global = Bounds(1e-10, 1e15); 	// in 1e48 ergs
	B.zBounds_global = Bounds(1e-3, 10);  B.zGridLen = 30;
	B.kBounds_global = Bounds(1e-2, 1e3); B.kGridLen = 200;
	B.SBounds_global = Bounds(1e-20, 1e-4);  B.SGridLen = 20;
	B.EBounds_global = Bounds(0.1_GeV, 900._GeV);  B.EGridLen = 40;
	B.GammaGridLen = 13;
	
	B.IntensityBins = {Bounds(0.1_GeV ,.14_GeV), Bounds(0.14_GeV ,0.2_GeV), Bounds(0.2_GeV ,0.28_GeV), Bounds(0.28_GeV ,0.4_GeV),
								Bounds(0.4_GeV ,0.57_GeV), Bounds(0.57_GeV ,0.8_GeV),  Bounds(0.8_GeV ,1.1_GeV), Bounds(1.1_GeV ,1.6_GeV), 
								Bounds(1.6_GeV, 2.3_GeV), Bounds(2.3_GeV, 3.2_GeV), Bounds(3.2_GeV, 4.5_GeV), Bounds(4.5_GeV, 6.4_GeV),
								Bounds(6.4_GeV, 9.1_GeV), Bounds(9.1_GeV, 13._GeV), Bounds(13._GeV, 18._GeV), Bounds(18._GeV, 26._GeV),
								Bounds(26._GeV, 36._GeV), Bounds(36._GeV, 51._GeV), Bounds(51._GeV, 72._GeV), Bounds(72._GeV, 100._GeV),
								Bounds(100._GeV, 140._GeV), Bounds(140._GeV, 200._GeV), Bounds(200._GeV, 290._GeV), Bounds(290._GeV, 410._GeV),
								Bounds(410._GeV, 580._GeV), Bounds(580._GeV, 820._GeV)};
	
	B.APSBins = {Bounds(0.5_GeV, 0.72_GeV ), Bounds(0.72_GeV, 1.04_GeV), Bounds(1.04_GeV, 1.38_GeV), Bounds(1.38_GeV, 1.99_GeV), 
					Bounds(1.99_GeV, 3.15_GeV), Bounds(3.15_GeV, 5.0_GeV), Bounds(5.0_GeV, 7.23_GeV), Bounds(7.23_GeV, 10.24_GeV), 
					Bounds(10.24_GeV, 21.83_GeV), Bounds(21.83_GeV, 50.0_GeV), Bounds(50.0_GeV, 95.27_GeV), Bounds(95.27_GeV, 199.05_GeV),
					Bounds(199.05_GeV, 500.0_GeV)  };

	std::vector<std::shared_ptr<AstrophysicalSource> > AstrophysicalSources;
	std::vector<std::shared_ptr<DarkMatter> > dmModels;
	
	
	//AstrophysicalSources.push_back(std::make_shared<MAGN>(CM, tau));
	//AstrophysicalSources.push_back(std::make_shared<FSRQ>(CM, tau));
	//AstrophysicalSources.push_back(std::make_shared<LISP>(CM, tau));
	//AstrophysicalSources.push_back(std::make_shared<HSP>(CM, tau));
	//AstrophysicalSources.push_back(std::make_shared<NGSFG>(CM, tau));
	//B.calculateIntensityAndAPSForAstrophysicalSources(AstrophysicalSources); 
	
	dmModels.push_back(std::make_shared<AnnihilatingDM>(CM, HM, tau, dNdLogx, 1e6, 3e-26));
	dmModels.push_back(std::make_shared<DecayingDM>(CM, HM, tau, dNdLogx, 10, 1.4e17));
	
	B.calculateIntensityAndAutocorrelationForDM(dmModels);
	
	//if(plot) B.SavePlots("plots/");
	//if(plot) plotIntensities(AstrophysicalSources, "plots/");
	
	
	
	
	for(unsigned int i = 0; i < dmModels.size(); i++)
	{
		std::cout << "DM Intensity: " ;
		for(unsigned int j = 0; j < dmModels.at(i)->Intensity.size(); j++)
			std::cout<< dmModels.at(i)->Intensity.at(j) << '\t';
		std::cout << std::endl;
	}

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

void plotIntensities(std::vector<std::shared_ptr<AstrophysicalSource> > sources, std::string path)
{
	Canvas IntensityCanvas("Int", "Intensity");
	IntensityCanvas().SetLogx(1); IntensityCanvas().SetLogy(1); 
	for(unsigned int i = 0; i < sources.size(); i++)
	{
		//TGraph* g = new TGraph(sources[i]->MakeGraph());
		//g->SetLineColor(i+1);
		//IntensityCanvas.AddGraph(g, sources[i]->Name, "L", true);
	}
	IntensityCanvas.Draw("A");
	IntensityCanvas().BuildLegend();
	IntensityCanvas().SaveAs((path + "ASIntensity.jpg").c_str());
}


