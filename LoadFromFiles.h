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

#include "TROOT.h"
#include "TFile.h"

#include "Constants.h"
#include "InterpolationWrapper.h"
#include "EBLAbsorbtionCoefficient.h"
#include "LinearMatterPowerSpectrum.h"
#include "AngularPowerSpectrum.h"

// Loads the data for the EBL Absorption coefficient
// and returns a 2D interpolation object, depending on energy and redshift
std::shared_ptr<EBLAbsorbtionCoefficient> LoadEBLAbsorbtionCoefficient(std::fstream& file)
{
	assert(file.is_open());
	// check numbers
	int EGridSize = 49; int zGridSize = 500;
	double* EGrid = new double[EGridSize];
	double* zGrid = new double[zGridSize];
	double** tau = new double*[EGridSize]; for(int i = 0; i < EGridSize; i++) tau[i] = new double[zGridSize];
	
	// Get Energy data:
	for(int i = 0; i < EGridSize; i++) {	file >> EGrid[i]; EGrid[i]*=1e3_GeV;	}
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

// Loads the linear matter power spectrum from a list of files
// returns 2D interpolation object depending on wavenumber and resdshift
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

// loads the energy spectrum of dark matter from simulations
// returns 2D interpolation object depending on mass and logX, where x = 2E/m 
// check papers
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

// Loads the energy bins and intensity measurements from a single file
void LoadDGRB(std::fstream& file, std::vector<Bounds>& IntensityBins, std::vector<Measurement>& DGRB)
{
	std::string line; std::stringstream ss; 
	Bounds b;  Measurement im;
	while(!file.eof())
	{
		std::getline(file, line);
		if(line.empty()) continue;
		ss = std::stringstream(line);
		ss >> b.first; ss >> b.second; IntensityBins.push_back(b);
		ss >> im.first; ss >> im.second; DGRB.push_back(im);
	}
}

// Loads the whole APS (Cl for different energy bins and multipole)
std::shared_ptr<AngularPowerSpectrum<Measurement> > LoadAnistropy(std::string directory, int nBin, int nMul, std::vector<double>& Multipoles)
{
	std::string line; std::stringstream ss; double buf; Measurement m;
	auto APS = std::make_shared<AngularPowerSpectrum<Measurement> >(nBin, nBin, nMul);
	for( int i = 0; i < nBin; i++)
	{
		for( int j = i; j < nBin; j++)
		{
			std::fstream file(directory + "Cl_p7_2FGL_default_" + std::to_string(i+1) + "_" + std::to_string(j+1) + ".txt", std::fstream::in);
			assert(file.is_open());
			std::getline(file, line);	// ignore first line
			int k = 0;
			while(!file.eof())
			{
				std::getline(file, line);
				if(line.empty()) continue;
				ss = std::stringstream(line);
				ss >> buf; ss >> buf; //ignore
				ss >> buf; if(i ==1 && j == 1) Multipoles.push_back(buf);  
				ss >> buf; m.first = buf; 	// value
				ss >> buf; m.second = buf;	// error
				//std::cout << i << ", " << j << ", " << k << ": " << m.first << " +- " << m.second << std::endl;
				APS->at(i, j, k) = m;
				if(i != j) APS->at(j, i, k) = m;	// symmetry
				k += 1;
			}
		}
	}
	return APS;
}

// Only loads Cp data
std::shared_ptr<AngularPowerSpectrum<Measurement> > LoadAnistropy(std::string filepath, int nBin)
{
	std::string line; std::stringstream ss;
	std::fstream file(filepath, std::fstream::in);
	assert(file.is_open());
	auto APS = std::make_shared<AngularPowerSpectrum<Measurement> >(nBin, nBin, 1);
	std::getline(file, line); // ignore first line
	for( int i = 0; i < nBin; i++)
	{
		std::getline(file, line);
		ss = std::stringstream(line);
		for( int j = 0; j < nBin; j++)
		{
			ss >> APS->at(i, j, 0).first;
		}
	}
	std::getline(file, line);
	for( int i = 0; i < nBin; i++)
	{
		std::getline(file, line);
		ss = std::stringstream(line);
		for( int j = 0; j < nBin; j++)
		{
			ss >> APS->at(i, j, 0).second;
		}
	}
	/*for(int i = 0; i < nBin; i++)
	{
		for( int j = 0; j <= i; j++)
			std::cout << i << ", " << j << ": " << APS->at(i, j, 0).first << " +- " <<  APS->at(i, j, 0).second << std::endl;
	}*/
	return APS;
}
#endif
