#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <tuple>

// Mass reference
constexpr double operator "" _M_solar(long double d) { return d; }

// Energy reference
const double GeV = 1;
constexpr  double operator "" _GeV(long double d) { return d; }

// Flux reference
constexpr double operator "" _photonspercm2s(long double d) { return d; }

// Luminosity reference
constexpr double operator "" _ergpers(long double d) { return d*1e-48; }

// Solar Luminosity
const double L_solar = 3.85e33_ergpers;

// Distance reference 
constexpr double operator "" _Mpc(long double d) { return d; }

// Speed of light km/s
const double c_0 = 2.9979e5;

// Hubble constant in km / s*MPc
const double H_0 = 67.74/1._Mpc;

// Hubble parameter
const double h = H_0/100.;



typedef std::pair<double, double> Bounds ;			// To store Integration Bounds (left, right)
typedef std::pair<double, double> Measurement;  // value and error

#endif
