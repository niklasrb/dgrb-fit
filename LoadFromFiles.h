#ifndef LOADFROMFILES_H
#define LOADFROMFILES_H
/// Header file for collecting functions that load data from a file and interpolate them
/// Mind the different file structures!

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <sstream>

#include "InterpolationWrapper.h"
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

std::shared_ptr<gsl2DInterpolationWrapper> LoaddNdLogx(std::fstream& file)
{
	std::string line;
	std::stringstream ss;
	double buf;
	std::vector<double> mass, logx, bbar;
	
	std::getline(file, line);	// first line is text
	
	while(!file.eof())
	{
		std::getline(file, line);
		if(line.empty()) continue;
		ss = std::stringstream(line);
		ss >> buf;	// first column is mass
		if(mass.size() == 0) mass.push_back(buf);
		else if(mass[mass.size()-1] != buf) mass.push_back(buf);
		
		ss >> buf; 	// second column is logx
		if(logx.size() == 0) logx.push_back(buf);
		else if(logx[logx.size()-1] < buf) logx.push_back(buf); // logx is in increasing order
		
		for(unsigned int i=0; i < 11; i++) ss >> buf;	// 11 uninteresting columns
		ss >> buf; // b case
		bbar.push_back(buf);		
	}

	auto NofLogx = std::make_shared<gsl2DInterpolationWrapper>(mass.data(), mass.size(), logx.data(), logx.size());
	for(unsigned int i = 0; i < mass.size(); i++)
	{
		for(unsigned int j = 0; j < logx.size(); j++)
				NofLogx->Val(i,j) = bbar.at(i*logx.size() + j);
	}
	NofLogx->Initialize();
	return NofLogx;
}

#endif
