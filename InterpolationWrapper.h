#ifndef INTERPOLATIONWRAPPER_H
#define INTERPOLATIONWRAPPER_H

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "TROOT.h"
#include "TGraph.h"
#include "Constants.h"
#include <cassert>
#include <iostream>
#include <istream>
#include <ostream>
#include <cmath>
#include <algorithm>
#include <memory>

class gsl1DInterpolationWrapper;
class gsl2DInterpolationWrapper;
class Interpolation3DWrapper;


/// 1D Interpolation
class gsl1DInterpolationWrapper
{
friend class gsl2DInterpolationWrapper;
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
	gsl1DInterpolationWrapper(std::istream& i, const gsl_interp_type* T, double outOfBounds);
	gsl1DInterpolationWrapper();
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
	
	// Saves data in the stream
	void Save(std::ostream& o);
	
	TGraph MakeGraph();
	
	// Allows multiplication if they have the same grid!
	friend gsl1DInterpolationWrapper operator*(const gsl1DInterpolationWrapper& iw1, const gsl1DInterpolationWrapper& iw2)
	{
		assert(iw1.n == iw2.n);
		for(unsigned int i = 0; i < iw1.n; i++) assert(iw1.x[i] == iw2.x[i]);
		double* y = new double[iw1.n];
		for(unsigned int i = 0; i < iw1.n; i++) y[i] = iw1.y[i] * iw2.y[i];
		gsl1DInterpolationWrapper res(iw1.x, iw1.n, y, iw1.T, iw1.outOfBounds*iw2.outOfBounds);
		delete []y;
		return res;
	}
	
	Bounds get_xBounds() { return xBounds; }
};


/// 2D Interpolation
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
	bool initiated;

public:
	// Main constructor
	gsl2DInterpolationWrapper(const double* _x, const unsigned int n_x, const double* _y, const unsigned int n_y, const double** _z, const gsl_interp2d_type* T);
	gsl2DInterpolationWrapper(const double* _x, const unsigned int n_x, const double* _y, const unsigned int n_y, const gsl_interp2d_type* T);
	gsl2DInterpolationWrapper(std::istream& in, const gsl_interp2d_type* T);
	
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
	
	// Allows acces to z
	double& Val(unsigned int i, unsigned int y);
	
	// Performs initiation of gsl if values weren't specified on construction
	void Initialize();
	
	// Outputs most data for debugging purposes
	void print();
	
	// Saves the data in the stream
	void Save(std::ostream& o);
	
	Bounds get_xBounds() { return xBounds; }
	Bounds get_yBounds() { return yBounds; }
	
	static gsl2DInterpolationWrapper Combine(const std::vector<std::shared_ptr<gsl1DInterpolationWrapper> > &splines, const std::vector<double>& y);
};


/// gsl1DInterpolation


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

// Constructs an Interpolator that always return NAN
gsl1DInterpolationWrapper::gsl1DInterpolationWrapper() : T(gsl_interp_linear), outOfBounds(NAN)
{
	n = 2;
	x = new double[n]; y = new double[n];
	x[0] = 0; x[1] = 1;  y[0] = NAN; y[1] = NAN;
	xBounds.first = 0; xBounds.second = -1;
	spline = gsl_interp_alloc(T, n);
	x_acc = gsl_interp_accel_alloc();
	gsl_interp_init(spline, x, y, n);
	
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
	std::cout << "gsl1D interp wrapper: " << spline;
	std::cout << '\t' << n <<" values: " << std::endl;
	for(unsigned int i = 0; i < n; i++)
			std::cout << "[" << x[i] <<  ", " << y[i] << "]" << (i == (n-1) ? '\n' : '\t');
}

void gsl1DInterpolationWrapper::Save(std::ostream& o)
{
	o << n << std::endl;
	for(unsigned int i = 0; i < n; i++) o << x[i] << '\t' << y[i] << std::endl;
}

gsl1DInterpolationWrapper::gsl1DInterpolationWrapper(std::istream& in, const gsl_interp_type* T= gsl_interp_linear, double outOfBounds = NAN) : T(T), outOfBounds(outOfBounds)
{
	in >> n; assert(n > 0);
	x = new double[n];  y = new double[n];
	spline = gsl_interp_alloc(T, n);
	x_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n; i++) in >> x[i] >> y[i]; 
	gsl_interp_init(spline, x, y, n);
	xBounds.first = x[0];  xBounds.second = x[n-1];	
}


TGraph gsl1DInterpolationWrapper::MakeGraph()
{
	return TGraph(n, x, y);
}


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
			if(i ==0) y[j] = _y[j];
			gsl_interp2d_set(spline, z, i, j, _z[i][j]);
			//z[i+n_x*j] = _z[i][j];
		}
	}
	gsl_interp2d_init(spline, x, y, z, n_x, n_y);
	xBounds.first = x[0];  xBounds.second = x[n_x-1];
	yBounds.first = y[0];  yBounds.second = y[n_y-1];
	initiated = true;
}

gsl2DInterpolationWrapper::gsl2DInterpolationWrapper(const double* _x, const unsigned int n_x, const double* _y, const unsigned int n_y, const gsl_interp2d_type* T= gsl_interp2d_bilinear) : n_x(n_x), n_y(n_y), T(T)
{
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_x*n_y];
	spline = gsl_interp2d_alloc(T, n_x, n_y);
	x_acc = gsl_interp_accel_alloc();
	y_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n_x; i++)
		x[i] = _x[i];
	for(unsigned int j = 0; j < n_y; j++)
		y[j] = _y[j];
	xBounds.first = x[0];  xBounds.second = x[n_x-1];
	yBounds.first = y[0];  yBounds.second = y[n_y-1];
	initiated = false;
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
	initiated = true;
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
	std::swap(T, tmp.T); std::swap(xBounds, tmp.xBounds); std::swap(yBounds, tmp.yBounds);  std::swap(tmp.initiated, initiated);
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
	std::cout << "gsl2D interp wrapper: " << spline << '\t';
	std::cout << n_x << "*" << n_y << " values: " << std::endl;
	for(unsigned int i = 0; i < n_x; i++)
		for(unsigned int j = 0; j < n_y; j++)
			std::cout << "[" << x[i] << ", " << y[j] << ", " << z[j*n_x + i] << "]" << (i*j == (n_x-1)*(n_y-1) ? '\n' : '\t');
}

void gsl2DInterpolationWrapper::Save(std::ostream& o)
{
	o << n_x << '\t' << n_y << std::endl;
	for(unsigned int i = 0; i < n_x; i++)	// first write x value in a line
		o << x[i] << (i < n_x-1 ? '\t' : '\n');
	for(unsigned int i = 0; i < n_y; i++)	// then y values
		o << y[i] << (i < n_y-1 ? '\t' : '\n');
	for(unsigned int i = 0; i < n_x; i++)	// then write z values in one line
	{
		for(unsigned int j = 0; j < n_y; j++)
		{
			o << z[i+n_x*j] << ((i < n_x-1)||(j < n_y-1) ? '\t' : '\n');
		}
	}
}

gsl2DInterpolationWrapper::gsl2DInterpolationWrapper(std::istream& in, const gsl_interp2d_type* T= gsl_interp2d_bilinear) : T(T)
{
	in >> n_x >> n_y;
	assert(n_x >= 2 && n_y >= 2);
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_x*n_y];
	spline = gsl_interp2d_alloc(T, n_x, n_y);
	x_acc = gsl_interp_accel_alloc();
	y_acc = gsl_interp_accel_alloc();
	for(unsigned int i = 0; i < n_x; i++)
		in >> x[i];
	for(unsigned int i = 0; i < n_y; i++)
		in >> y[i];
	for(unsigned int i=0; i < n_x; i++)
	{
		for(unsigned int j = 0; j < n_y; j++)
		{
			in >> z[i+n_x*j];
		}
	}
	gsl_interp2d_init(spline, x, y, z, n_x, n_y);
	xBounds.first = x[0];  xBounds.second = x[n_x-1];
	yBounds.first = y[0];  yBounds.second = y[n_y-1];
	initiated = true;
}

double& gsl2DInterpolationWrapper::Val(unsigned int i, unsigned int j)
{
	assert(i < n_x && j < n_y);
	return z[i+n_x*j];
}
	
	// Performs initiation of gsl
void gsl2DInterpolationWrapper::Initialize()
{
	if(initiated) return;
	gsl_interp2d_init(spline, x, y, z, n_x, n_y);
	initiated = true;
}

gsl2DInterpolationWrapper gsl2DInterpolationWrapper::Combine(const std::vector<std::shared_ptr<gsl1DInterpolationWrapper> >& splines, const std::vector<double>& y)
{
	assert(splines.size() >= 2);
	assert(y.size() == splines.size());
	unsigned int n_y = y.size();
	unsigned int n_x = splines[0]->n;
	
	for(unsigned int i = 1; i < splines.size(); i++)
	{
		assert(splines[i]->n == splines[i-1]->n);
		for(unsigned int j = 0; j < splines[i]->n; j++) assert(splines[i]->x[j] == splines[i-1]->x[j]);
	}
	gsl2DInterpolationWrapper iw(splines[0]->x, n_x, y.data(), n_y);
	for(unsigned int i = 0; i < n_x; i++)
	{
		for(unsigned int j = 0; j < n_y; j++)
		{
			iw.Val(i, j) = splines[j]->y[i]; // remember that y of gsl1DIW is the third axis in this case
		}
	}
	iw.Initialize();
	return iw;
}


/// 3D interpolation

class Interpolation3DWrapper
{
protected:
double* x; unsigned int n_x; 
double* y; unsigned int n_y; 
double* z; unsigned int n_z; 
double* f;

double interp(const double& _x, const double& _y, const double& _z);	// no bounds checking!
double extrap(const double& _x, const double& _y, const double& _z);
double& access(unsigned int i, unsigned int j, unsigned int k);			// no bounds checking!

public:

Interpolation3DWrapper(const double* _x,unsigned int n_x, const double* _y, unsigned int n_y, const double* z,unsigned int n_z);
Interpolation3DWrapper(std::istream& in);
Interpolation3DWrapper(const Interpolation3DWrapper& iw);
~Interpolation3DWrapper();
Interpolation3DWrapper& operator =(const Interpolation3DWrapper& iw);
double& Val(unsigned int i, unsigned int j, unsigned int k);			// with bounds checking
double Eval(const double& _x, const double& _y, const double& _z);		// with bounds checking
void print();
void Save(std::ostream& o);

};

Interpolation3DWrapper::Interpolation3DWrapper(const double* _x,unsigned int n_x, const double* _y,unsigned int n_y, const double* _z,unsigned int n_z) : n_x(n_x), n_y(n_y), n_z(n_z)
{
	assert(n_x > 1 && n_y > 1 && n_z > 1);
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_z];
	
	for(unsigned int i = 0; i < n_x; i++)	
	{
		x[i] = _x[i];
		if(i > 0) assert(x[i] > x[i-1]);
	}
	for(unsigned int i = 0; i < n_y; i++)
	{
		y[i] = _y[i];
		if(i > 0) assert(y[i] > y[i-1]);
	}
	for(unsigned int i = 0; i < n_z; i++)
	{
		z[i] = _z[i];
		if(i > 0) assert(z[i] > z[i-1]);
	}
	
	f = new double[n_x*n_y*n_z];
	for(unsigned int i = 0; i < n_x*n_y*n_z; i++) f[i] = 0;
}

Interpolation3DWrapper::~Interpolation3DWrapper()
{
	delete []x;
	delete []y;
	delete []z;
	delete []f;
}

double& Interpolation3DWrapper::Val(unsigned int i,unsigned int j,unsigned int k)
{
	assert(i < n_x && j < n_y && k < n_z);
	return access(i, j, k);
}

double Interpolation3DWrapper::interp(const double& _x, const double& _y, const double& _z)
{
	unsigned int i = 1; unsigned int j = 1; unsigned int k = 1;
	while(x[i] < _x && i <= n_x-2) i++;
	while(y[j] < _y && j <= n_y-2) j++;
	while(z[k] < _z && k <= n_z-2) k++;
	
	const double x_d = (_x - x[i-1])/(x[i] - x[i-1]);
	const double y_d = (_y - y[j-1])/(y[j] - y[j-1]);
	const double z_d = (_z - z[k-1])/(z[k] - z[k-1]);
	
	// interpolate in x
	const double c00 = access(i-1, j-1, k-1)*(1.-x_d) + access(i, j-1, k-1)*x_d;
	const double c01 = access(i-1, j-1, k)*(1.-x_d) + access(i, j-1, k)*x_d;
	const double c10 = access(i-1, j, k-1)*(1.-x_d) + access(i, j, k-1)*x_d;
	const double c11 = access(i-1, j, k)*(1.-x_d) + access(i, j, k)*x_d;
	
	// interpolate in y
	const double c0 = c00*(1.-y_d) + c10*y_d;
	const double c1 = c01*(1.-y_d) + c11*y_d;
	
	return c0*(1.-z_d) + c1*z_d;
}

double Interpolation3DWrapper::Eval(const double& _x, const double& _y, const double& _z)
{
	assert(x[0] <= _x && x[n_x-1] >= _x);
	assert(y[0] <= _y && y[n_y-1] >= _y);
	assert(z[0] <= _z && z[n_z-1] >= _z);
	return interp(_x, _y, _z);
}

double& Interpolation3DWrapper::access(unsigned int i, unsigned int j, unsigned int k)
{
	return f[i + n_x * (j + n_y * k)];	// god please be right
}

Interpolation3DWrapper::Interpolation3DWrapper(const Interpolation3DWrapper& iw) : n_x(iw.n_x), n_y(iw.n_y), n_z(iw.n_z)
{
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_z];
	f = new double[n_x*n_y*n_z];
	for(unsigned int i = 0; i < n_x; i++) x[i] = iw.x[i];
	for(unsigned int i = 0; i < n_y; i++) y[i] = iw.y[i];
	for(unsigned int i = 0; i < n_z; i++) z[i] = iw.z[i];
	for(unsigned int i = 0; i < n_x*n_y*n_z; i++) f[i] = iw.f[i];
}

Interpolation3DWrapper& Interpolation3DWrapper::operator =(const Interpolation3DWrapper& iw)
{
	if(this == &iw) return *this;
	Interpolation3DWrapper tmp(iw);
	std::swap(x, tmp.x); std::swap(n_x, tmp.n_x);
	std::swap(y, tmp.y); std::swap(n_y, tmp.n_y);
	std::swap(z, tmp.z); std::swap(n_z, tmp.n_z);
	std::swap(f, tmp.f);	
	return *this;
}

void Interpolation3DWrapper::print()
{
	std::cout << "3D interp " << n_x << '*' << n_y << '*' << n_z << " values: ";
	for(unsigned int i = 0; i < n_x; i++)
	{
		for(unsigned int j = 0; j < n_y; j++)
		{
			for(unsigned int k = 0; k < n_z; k++)
			{
				std::cout << '[' << x[i] << ", " << y[j] << ", " << z[k] << ", " << access(i, j, k) << "]" << '\t';
			}
		}
	}
	std::cout << std::endl;
}

void Interpolation3DWrapper::Save(std::ostream& o)
{
	o << n_x << '\t' << n_y << '\t' << n_z << '\n';
	for(unsigned int i = 0; i < n_x; i++)
		o << x[i] << (i < n_x-1 ? '\t' : '\n'); 
	for(unsigned int i = 0; i < n_y; i++)
		o << y[i] << (i < n_y-1 ? '\t' : '\n'); 
	for(unsigned int i = 0; i < n_z; i++)
		o << z[i] << (i < n_z-1 ? '\t' : '\n');
	for(unsigned int i = 0; i < n_x; i++)
	{
		for(unsigned int j = 0; j < n_y; j++)
		{
			for(unsigned int k = 0; k < n_z; k++)
			{
				o << access(i,j,k) << ( (i<n_x-1)||(j<n_y-1)||(k<n_z-1) ? '\t' : '\n');
			}
		}
	} 
}
Interpolation3DWrapper::Interpolation3DWrapper(std::istream& in) 
{
	in >> n_x >> n_y >> n_z;
	assert(n_x > 1 && n_y > 1 && n_z > 1);
	x = new double[n_x];
	y = new double[n_y];
	z = new double[n_z];
	f = new double[n_x*n_y*n_z];
	for(unsigned int i = 0; i < n_x; i++)
		in >> x[i];
	for(unsigned int i = 0; i < n_y; i++)
		in >> y[i];
	for(unsigned int i = 0; i < n_z; i++)
		in >> z[i];
	for(unsigned int i = 0; i < n_x; i++)
	{
		for(unsigned int j = 0; j < n_y; j++)
		{
			for(unsigned int k = 0; k < n_z; k++)
			{
				in >> access(i,j,k);
			}
		}
	} 
}



#endif
