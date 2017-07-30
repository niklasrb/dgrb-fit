#ifndef ANGULARPOWERSPECTRUM_H
#define ANGULARPOWERSPECTRUM_H

#include <cmath>
#include <vector>
#include <cassert>

#include "Constants.h"

template<typename T> 
class AngularPowerSpectrum
{
protected:
	T*** data;		// 3dim Matrix to hold the APS
	unsigned int nBin1, nBin2, nMul;	// dimension of this matrix
	//std::vector<Bounds> EBins1, EBins2;
	//std::vector<int> Multipole;
	
public:
	virtual T& operator ()(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)
	{
		assert(EBin1 < nBin1  && EBin2 < nBin2 && Multipole < nMul);
		return data[EBin1][EBin2][Multipole];
	}
	
	/// Use this constructor if you dont't have the data yet, or don't want to deal with memory
	AngularPowerSpectrum(unsigned int nBin1, unsigned int nBin2, unsigned int nMul) : nBin1(nBin1), nBin2(nBin2), nMul(nMul)
	{
		data = new T**[nBin1];
		for(unsigned int i = 0; i < nBin1; i++)
		{
			data[i] = new T*[nBin2];
			for(unsigned int j = 0; j < nBin2; j++)
			{
				data[i][j] = new T[nMul];
			}
		}
	}
	
	/// Use this constructor if you already have the data in the correct format
	AngularPowerSpectrum(unsigned int nBin1, unsigned int nBin2, unsigned int nMul, T*** APS) : AngularPowerSpectrum(nBin1, nBin2, nMul)
	{
		for(unsigned int i = 0; i < nBin1; i++)
		{
			for(unsigned int j = 0; j < nBin2; j++)
			{
				for(unsigned int k = 0; k < nMul; k++)
				{
					data[i][j][k] = APS[i][j][k];
				}
			}
		}
	}
	
	~AngularPowerSpectrum()
	{
		for(unsigned int i = 0; i < nBin1; i++)
		{
			for(unsigned int j = 0; j < nBin2; j++)
			{
				delete []data[i][j];
			}
			delete []data[i];
		}
		delete []data;
	}
	
	AngularPowerSpectrum(const AngularPowerSpectrum& aps) = delete;				// forbid copying for now
	AngularPowerSpectrum operator =(const AngularPowerSpectrum& aps) = delete;
	
	
};


/// This class makes use of the fact that C_p is constant in Multipole
template<typename T> 
class AstrophysicalSourceAPS :	public AngularPowerSpectrum<T>
{
	AstrophysicalSourceAPS(std::vector<Bounds>& EBins) : AngularPowerSpectrum<T>::AngularPowerSpectrum(EBins.size(), EBins.size(), 1)	{	}
	
	
	T& operator ()(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole) override
	{
		assert(EBin1 < this->nBin1  && EBin2 < this->nBin2);
		return this->data[EBin1][EBin2][0];
	}
	
};
/*
/// Specific implementation that is useful when you have the APS as a list of doubles
template<> 
class AstrophysicalSourceAPS<double>
{
	AstrophysicalSourceAPS(std::vector<Bounds>& EBins, std::vector<double>& C_p) : AstrophysicalSourceAPS<double>::AstrophysicalSourceAPS(EBins)
	{
		assert(EBins.size() == C_p.size());
		
		for(unsigned int i = 0; i < EBin.size(); i++)
		{
			for(unsigned int j = 0; j < EBin.size(); j++)
			{
				this->data[i][j][0] = sqrt(C_p[i]*C_p[j]);
			}
		}
	}
};
*/


#endif
