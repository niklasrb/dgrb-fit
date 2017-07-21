#ifndef LOADFROMFILES_H
#define LOADFROMFILES_H
/// Header file for collecting functions that load data from a file and interpolate them
/// Mind the different file structures!

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>

#include "gsl2DInterpolationWrapper.h"
#include "EBLAbsorbtionCoefficient.h"
#include "LinearMatterPowerSpectrum.h"

std::shared_ptr<EBLAbsorbtionCoefficient> LoadEBLAbsorbtionCoefficient(std::fstream& file)
{
	assert(file.is_open());
	// check numbers
	int EGridSize = 49; int zGridSize = 500;
	double* EGrid = new double[EGridSize];
	double* zGrid = new double[zGridSize];
	double** tau = new double*[EGridSize]; for(int i = 0; i < EGridSize; i++) tau[i] = new double[zGridSize];
	
	// Get Energy data:
	for(int i = 0; i < EGridSize; i++) file >> EGrid[i];
	// Read per Block
	for(int j = 0; j < zGridSize; j++)
	{
		file >> zGrid[j];
		//std::cout << j << ": reading redshift " << zGrid[j] << std::endl;
		for(int i = 0; i < EGridSize; i++)
			file >> tau[i][j];
	}
	gsl2DInterpolationWrapper spline(EGrid, EGridSize, zGrid, zGridSize, (const double**) tau);
	delete []EGrid; delete []zGrid; for(int i = 0; i < EGridSize; i++) delete []tau[i]; delete []tau;
	return std::make_shared<EBLAbsorbtionCoefficient>(spline);
}


std::shared_ptr<LinearMatterPowerSpectrum> LoadLinearMatterPowerSpectrum(std::vector<std::fstream>& files)
{
	unsigned int kGridSize = 643;
	unsigned int zGridSize = files.size();
	
	double* z = new double[zGridSize];
	double* k = new double[kGridSize];
	double** P = new double*[kGridSize]; for(unsigned int i = 0; i < kGridSize; i++) P[i] = new double[zGridSize]; 
	
	for(unsigned int j = 0; j < files.size(); j++)
	{
		files.at(j) >> z[j];
		for(unsigned int i = 0; i < kGridSize; i++)
		{
			files.at(j) >> k[i];
			files.at(j) >> P[i][j];
		}
	}
	
	gsl2DInterpolationWrapper spline(k, kGridSize, z, zGridSize, (const double**)P);
	delete []z; delete []k;
	for(unsigned int i = 0; i < kGridSize; i++) delete []P[i];
	delete []P;
	//spline.print();
	return std::make_shared<LinearMatterPowerSpectrum>(spline);
}



#endif
