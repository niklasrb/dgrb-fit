#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <tuple>

// Speed of light km/s
const double c_0 = 2.9979e5;

// Hubble constant in km / s*MPc
const double H_0 = 67.74;

// Hubble parameter
const double h = H_0/100.;

// Gravtity constant in km^3 / kg*s^2
const double G = 6.67e-20;

// Solar Mass in kg
const double M_solar = 1.988e20;

// Joule To GeV
const double e = 1.6e-13;


// Mass reference
//constexpr double operator "" _M_solar(long double d) {return d;}

// Energy reference
const double GeV = 1;
constexpr  double operator "" _GeV(long double d) { return d; }
constexpr  double operator "" _GeV(unsigned long long d) { return d; }

// Flux reference
constexpr double operator "" _photonspercm2s1(long double d) { return d; }

// Luminosity reference
constexpr double operator "" _ergpers(long double d) { return d; }

typedef std::pair<double, double> Bounds ;			// To store Integration Bounds (left, right)




#endif
