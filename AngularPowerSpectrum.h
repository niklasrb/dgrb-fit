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
	T* data;		// 3dim Matrix to hold the APS
	unsigned int nBin1, nBin2, nMul;	// dimension of this matrix
	
	T& access(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)
	{
		return data[EBin1 + nBin1 * (EBin2 + nBin2 * Multipole)];
	}
	
public:
	virtual T& operator ()(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)
	{
		assert(EBin1 < nBin1  && EBin2 < nBin2 && Multipole < nMul);
		return access(EBin1, EBin2, nMul);
	}
	
	unsigned int Bin1Size() { return nBin1; }
	unsigned int Bin2Size() { return nBin2; }
	unsigned int MultipoleNumber() { return nMul; }
	
	
	virtual T& at(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)
	{
		return (*this)(EBin1, EBin2, Multipole);
	}
	
	/// Use this constructor if you dont't have the data yet, or don't want to deal with memory
	AngularPowerSpectrum(unsigned int nBin1, unsigned int nBin2, unsigned int nMul) : nBin1(nBin1), nBin2(nBin2), nMul(nMul)
	{
		data = new T[nBin1*nBin2*nMul];
	}
	
	~AngularPowerSpectrum()
	{
		delete []data;
	}
	
	AngularPowerSpectrum(const AngularPowerSpectrum& aps) : nBin1(aps.nBin1), nBin2(aps.nBin2), nMul(aps.nMul)
	{
		data = new T[nBin1*nBin2*nMul];
		for(unsigned int i = 0; i < nBin1*nBin2*nMul; i++) data[i] = aps.data[i];
	}
	
	AngularPowerSpectrum& operator =(const AngularPowerSpectrum& aps)
	{
		if(this == &aps) return *this;
		AngularPowerSpectrum<T> tmp(aps);
		std::swap(nBin1, tmp.nBin1); std::swap(nBin2, tmp.nBin2);
		std::swap(nMul, tmp.nMul); std::swap(data, tmp.data);
		return *this;
	}
	
	
};


/// This class makes use of the fact that C_p is constant in Multipole
template<typename T> 
class AstrophysicalSourceAPS :	public AngularPowerSpectrum<T>
{
public:
	AstrophysicalSourceAPS(const std::vector<Bounds>& EBins) : AngularPowerSpectrum<T>::AngularPowerSpectrum(EBins.size(), EBins.size(), 1)	{	}
	
	
	T& operator ()(unsigned int EBin1, unsigned int EBin2)
	{
		assert(EBin1 < this->nBin1  && EBin2 < this->nBin2);
		return this->access(EBin1, EBin2, 0);
	}
	
	T& operator ()(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole) override
	{
		assert(EBin1 < this->nBin1  && EBin2 < this->nBin2);
		return this->access(EBin1, EBin2, 0);
	}
	
	T& at(unsigned int EBin1, unsigned int EBin2)
	{
		return (*this)(EBin1, EBin2);
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
