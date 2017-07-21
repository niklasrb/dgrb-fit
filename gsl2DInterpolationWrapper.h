#ifndef GSLINTERPOLATIONWRAPPER_H
#define GSLINTERPOLATIONWRAPPER_H

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Constants.h"
#include <iostream>
#include <cmath>
#include <algorithm>



class gsl2DInterpolationWrapper
{
protected:
	gsl_interp2d* spline;
	gsl_interp_accel* x_acc;
	gsl_interp_accel* y_acc;
	Bounds xBounds;
	Bounds yBounds;
	double* x; unsigned int n_x;
	double* y; unsigned int n_y;
	double* z;
	const gsl_interp2d_type* T;

public:
	// Main constructor
	gsl2DInterpolationWrapper(const double* _x, const unsigned int n_x, const double* _y, const unsigned int n_y, const double** _z, const gsl_interp2d_type* T);
	gsl2DInterpolationWrapper(const gsl2DInterpolationWrapper& w);
	~gsl2DInterpolationWrapper();
	
	gsl2DInterpolationWrapper& operator =(const gsl2DInterpolationWrapper& w);
	
	// calls Eval
	double operator ()(const double _x, const double _y);
	
	// Perfoms bounds checking and extrapolates if needed
	double Eval(const double _x, const double _y);
	
	// Perfoms bounds checking and returns false for out of bounds
	bool Interpolate(const double _x, const double _y, double* _z);
	
	// Perfoms bounds checking and returns NAN for out of bounds
	double Interpolate(const double _x, const double _y);
	
	// Extrapolates explicitly
	double Extrapolate(const double _x, const double _y);
	
	// Outputs most data for debugging purposes
	void print();
	
	Bounds get_xBounds() { return xBounds; }
	Bounds get_yBounds() { return yBounds; }
};




/// gsl2DInterpolation

gsl2DInterpolationWrapper::gsl2DInterpolationWrapper(const double* _x, const unsigned int n_x, const double* _y, const unsigned int n_y, const double** _z, const gsl_interp2d_type* T= gsl_interp2d_bilinear) : n_x(n_x), n_y(n_y), T(T)
{
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_x*n_y];
	spline = gsl_interp2d_alloc(T, n_x, n_y);
	x_acc = gsl_interp_accel_alloc();
	y_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n_x; i++)
	{
		x[i] = _x[i];
		for(unsigned int j = 0; j < n_y; j++)
		{
			y[j] = _y[j];
			gsl_interp2d_set(spline, z, i, j, _z[i][j]);
		}
	}
	gsl_interp2d_init(spline, x, y, z, n_x, n_y);
	xBounds.first = x[0];  xBounds.second = x[n_x-1];
	yBounds.first = y[0];  yBounds.second = y[n_y-1];
	
}

gsl2DInterpolationWrapper::gsl2DInterpolationWrapper(const gsl2DInterpolationWrapper& w) : xBounds(w.xBounds), yBounds(w.yBounds), n_x(w.n_x), n_y(w.n_y),  T(w.T)
{
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_x*n_y];
	spline = gsl_interp2d_alloc(T, n_x, n_y);
	x_acc = gsl_interp_accel_alloc();
	y_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n_x; i++)
	{
		x[i] = w.x[i];
		for(unsigned int j = 0; j < n_y; j++)
		{
			y[j] = w.y[j];
			gsl_interp2d_set(spline, z, i, j, w.z[i + n_x*j]);	// this is the way the values are stored internally, check documentation
		}
	}
	gsl_interp2d_init(spline, x, y, z, n_x, n_y);
}

gsl2DInterpolationWrapper::~gsl2DInterpolationWrapper()
{
	gsl_interp2d_free(spline);
	gsl_interp_accel_free(x_acc);
	gsl_interp_accel_free(y_acc);
	delete []x;
	delete []y;
	delete []z;
}

gsl2DInterpolationWrapper& gsl2DInterpolationWrapper::operator =(const gsl2DInterpolationWrapper& w)
{
	if(this == &w) return *this;
	gsl2DInterpolationWrapper tmp(w);
	std::swap(x, tmp.x); std::swap(n_x, tmp.n_x);
	std::swap(y, tmp.y); std::swap(n_y, tmp.n_y);
	std::swap(z, tmp.z);
	std::swap(spline, tmp.spline); std::swap(x_acc, tmp.x_acc); std::swap(y_acc, tmp.y_acc);
	std::swap(T, tmp.T); std::swap(xBounds, tmp.xBounds); std::swap(yBounds, tmp.yBounds);
	return *this;
}

double gsl2DInterpolationWrapper::operator ()(const double _x, const double _y)
{
	return Eval(_x, _y);	
}

double gsl2DInterpolationWrapper::Eval(const double _x, const double _y)
{
	if( _x < xBounds.first ||  _x > xBounds.second || _y < yBounds.first ||  _y > yBounds.second)  
		return gsl_interp2d_eval_extrap(spline, x, y, z, _x, _y, x_acc, y_acc);	
	return gsl_interp2d_eval(spline, x, y, z, _x, _y, x_acc, y_acc);	
}

bool gsl2DInterpolationWrapper::Interpolate(const double _x, const double _y, double* _z)
{
	if( _x < xBounds.first ||  _x > xBounds.second)  return false;
	if( _y < yBounds.first ||  _y > yBounds.second)  return false;
	return (EDOM != gsl_interp2d_eval_e(spline, x, y, z, _x, _y, x_acc, y_acc, _z));
}

double gsl2DInterpolationWrapper::Interpolate(const double _x, const double _y)
{
	if( _x < xBounds.first ||  _x > xBounds.second)  return NAN;
	if( _y < yBounds.first ||  _y > yBounds.second)  return NAN;
	return  gsl_interp2d_eval(spline, x, y, z, _x, _y, x_acc, y_acc);
}

double gsl2DInterpolationWrapper::Extrapolate(const double _x, const double _y)
{
	return gsl_interp2d_eval_extrap(spline, x, y, z, _x, _y, x_acc, y_acc);
}

void gsl2DInterpolationWrapper::print()
{
	std::cout << "gsl2D interp wrapper: " << spline << std::endl;
	std::cout << n_x << "*" << n_y << " values: " << std::endl;
	for(unsigned int i = 0; i < n_x; i++)
		for(unsigned int j = 0; j < n_y; j++)
			std::cout << "[" << x[i] << ", " << y[j] << ", " << z[j*n_x + i] << "]" << (i*j == (n_x-1)*(n_y-1) ? '\n' : '\t');
}


/// gsl1DInterpolation


class gsl1DInterpolationWrapper
{
protected:
	gsl_interp* spline;
	gsl_interp_accel* x_acc;
	Bounds xBounds;
	double* x; unsigned int n;
	double* y;
	const gsl_interp_type* T;
	const double outOfBounds;

public:
	// Main constructor
	gsl1DInterpolationWrapper(const double* _x, const unsigned int n, const double* _y, const gsl_interp_type* T, double outOfBounds);
	gsl1DInterpolationWrapper(const gsl1DInterpolationWrapper& w);
	~gsl1DInterpolationWrapper();
	
	gsl1DInterpolationWrapper& operator =(const gsl1DInterpolationWrapper& w);
	
	// calls Eval
	double operator ()(const double _x);
	
	// Perfoms bounds checking and return outOfBounds value in case
	double Eval(const double _x);
	
	// Perfoms bounds checking and returns false for out of bounds
	bool Interpolate(const double _x, double* _y);
	
	// Perfoms bounds checking and returns NAN for out of bounds
	double Interpolate(const double _x);
	
	// Performs bound checking (NAN for out of bounds) and returns the derivative at _x
	double Derivative(const double _x);
	
	// Outputs most data for debugging purposes
	void print();
	
	Bounds get_xBounds() { return xBounds; }
};


gsl1DInterpolationWrapper::gsl1DInterpolationWrapper(const double* _x, const unsigned int n, const double* _y, const gsl_interp_type* T= gsl_interp_linear, double outOfBounds = NAN) : n(n), T(T), outOfBounds(outOfBounds)
{
	x = new double[n];
	y = new double[n];
	spline = gsl_interp_alloc(T, n);
	x_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n; i++)
	{
		x[i] = _x[i];
		y[i] = _y[i];		
	}
	gsl_interp_init(spline, x, y, n);
	xBounds.first = x[0];  xBounds.second = x[n-1];	
}

gsl1DInterpolationWrapper::gsl1DInterpolationWrapper(const gsl1DInterpolationWrapper& w) : xBounds(w.xBounds), n(w.n),  T(w.T), outOfBounds(w.outOfBounds)
{
	x = new double[n];
	y = new double[n];
	spline = gsl_interp_alloc(T, n);
	x_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n; i++)
	{
		x[i] = w.x[i];
		y[i] = w.y[i];
	}
	gsl_interp_init(spline, x, y, n);
}

gsl1DInterpolationWrapper::~gsl1DInterpolationWrapper()
{
	gsl_interp_free(spline);
	gsl_interp_accel_free(x_acc);
	delete []x;
	delete []y;
}

gsl1DInterpolationWrapper& gsl1DInterpolationWrapper::operator =(const gsl1DInterpolationWrapper& w)
{
	if(this == &w) return *this;
	gsl1DInterpolationWrapper tmp(w);
	std::swap(x, tmp.x); std::swap(n, tmp.n);
	std::swap(y, tmp.y); 
	std::swap(spline, tmp.spline); std::swap(x_acc, tmp.x_acc);
	std::swap(T, tmp.T); std::swap(xBounds, tmp.xBounds);
	return *this;
}

double gsl1DInterpolationWrapper::operator ()(const double _x)
{
	return Eval(_x);	
}

double gsl1DInterpolationWrapper::Eval(const double _x)
{
	if( _x < xBounds.first ||  _x > xBounds.second)  
		return outOfBounds;	
	return gsl_interp_eval(spline, x, y, _x, x_acc);	
}

bool gsl1DInterpolationWrapper::Interpolate(const double _x, double* _y)
{
	if( _x < xBounds.first ||  _x > xBounds.second)  return false;
	return (EDOM != gsl_interp_eval_e(spline, x, y, _x, x_acc, _y));
}

double gsl1DInterpolationWrapper::Interpolate(const double _x)
{
	if( _x < xBounds.first ||  _x > xBounds.second)  return NAN;
	return  gsl_interp_eval(spline, x, y, _x, x_acc);
}

double gsl1DInterpolationWrapper::Derivative(const double _x)
{
	if( _x < xBounds.first ||  _x > xBounds.second)  return NAN;
	return gsl_interp_eval_deriv(spline, x, y, _x, x_acc);
}

void gsl1DInterpolationWrapper::print()
{
	std::cout << "gsl1D interp wrapper: " << spline << std::endl;
	std::cout << n <<" values: " << std::endl;
	for(unsigned int i = 0; i < n; i++)
			std::cout << "[" << x[i] <<  ", " << y[i] << "]" << (i == (n-1) ? '\n' : '\t');
}

#endif
