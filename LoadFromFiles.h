#ifndef LOADFROMFILES_H
#define LOADFROMFILES_H
/// Header file for collecting functions that load data from a file and interpolate them
/// Mind the different file structures!

#include <memory>
#include <fstream>
#include <string>
#include <cassert>

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


std::shared_ptr<LinearMatterPowerSpectrum> LoadLinearMatterPowerSpectrum(std::string file)
{
	double x[] = {0};
	double y[] = {0};
	double** z = new double*; *z = new double; **z = 0;
	gsl2DInterpolationWrapper spline(x, 1, y, 1, (const double**)z);
	delete *z; delete z;
	return std::make_shared<LinearMatterPowerSpectrum>(spline);
}



#endif
