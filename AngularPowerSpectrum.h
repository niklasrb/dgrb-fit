#ifndef ANGULARPOWERSPECTRUM_H
#define ANGULARPOWERSPECTRUM_H

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

#include "Constants.h"

/* Template class for holding different types of APS data
 * Models a 3D matrix with 2 energy bins and 1 multipole axis
 * The memory is constructed on constructor call and needs to be filled with () or at()
 */

template<typename T> 
class AngularPowerSpectrum
{
protected:
	T* data;		// 3dim Matrix to hold the APS
	unsigned int nBin1, nBin2, nMul;	// dimension of this matrix
	
	T& access(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)	// internal & without Bounds checking
	{
		return data[EBin1 + nBin1 * (EBin2 + nBin2 * Multipole)];
	}
	
public:
	virtual T& operator ()(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)	// Bounds checking
	{
		assert(EBin1 < nBin1  && EBin2 < nBin2 && Multipole < nMul);
		return access(EBin1, EBin2, Multipole);
	}
	
	unsigned int Bin1Size() { return nBin1; }
	unsigned int Bin2Size() { return nBin2; }
	unsigned int MultipoleNumber() { return nMul; }
	
	
	virtual T& at(unsigned int EBin1, unsigned int EBin2, unsigned int Multipole)
	{
		return (*this)(EBin1, EBin2, Multipole);
	}
	
	AngularPowerSpectrum(unsigned int nBin1, unsigned int nBin2, unsigned int nMul) : nBin1(nBin1), nBin2(nBin2), nMul(nMul)
	{
		data = new T[nBin1*nBin2*nMul];
	}
	
	~AngularPowerSpectrum()
	{
		delete [] data;
	}
	
	AngularPowerSpectrum(const AngularPowerSpectrum& aps) : nBin1(aps.nBin1), nBin2(aps.nBin2), nMul(aps.nMul)	// copy constructor
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
/// and makes essentially a 2D matrix
template<typename T> 
class AstrophysicalSourceAPS :	public AngularPowerSpectrum<T>
{
public:
	AstrophysicalSourceAPS(const std::vector<Bounds>& EBins) : AngularPowerSpectrum<T>::AngularPowerSpectrum(EBins.size(), EBins.size(), 1)	{	}
	
	AstrophysicalSourceAPS(unsigned int EBin1, unsigned int EBin2) : AngularPowerSpectrum<T>::AngularPowerSpectrum(EBin1, EBin2, 1)	{	}
	
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



#endif
