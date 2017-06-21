#ifndef SFG_H
#define SFG_H

#include "AstrophysicalSource.h"
#include "Constants.h"
#include <cmath>
#include <algorithm>

class SFG : public AstrophysicalSource
{
protected:

	
public:
	SFG(std::shared_ptr<CosmologyModel> _CM) : AstrophysicalSource(_CM, std::string("SFG"))
	{
		m_zBounds.first = 1; m_zBounds.second = 6;
		m_GammaBounds.first = ; m_GammaBounds.second = ;
	}
	
	double kCorrection(const double E, const double z, const double Gamma) override
	{

	}
	
	double EnergySpectrum(const double E, const double z, const double Gamma) override
	{
	}
	
	double GammaDistribution(const double Gamma) override
	{

	}			
	
	double LuminosityFunction(const double Luminosity, const double z, const double Gamma) override
	{

	}
	

};



#endif
