#include "Benchmarking.h"


int main(int argc, char** argv)
{
	auto CM = std::make_shared<LambdaCDM>();
	auto D = std::make_shared<FermiLAT>();
	Benchmark B(CM, D, true, true);
	//auto magn = std::make_shared<MAGN>(CM);
	//std::cout << magn->RescaledLuminosityFunction(1e9, 1, 2.37) << std::endl;
	
	B.AstrophysicalSources.push_back(std::make_shared<MAGN>(CM));
	B.EBins = {Bin(0.1_GeV, 0.5_GeV ,1._GeV), Bin(1._GeV, 5._GeV, 10._GeV), Bin(10._GeV, 50._GeV, 100._GeV), Bin(100._GeV, 500._GeV, 1000._GeV)};
	B.calculateIntensityAndAutocorrelation(5, 5);
	return 0;
}
