#ifndef CONSTANTS_H
#define CONSTANTS_H

//Speed of light km/s
const double c_0 = 2.9979e5;

//Hubble constant  km / s*MPc
const double H_0 = 67.74;

// Energy reference
const double GeV = 1;
constexpr  double operator "" _GeV(long double d) { return d; }
constexpr  double operator "" _GeV(unsigned long long d) { return d; }

// Flux reference
constexpr double operator "" _photonspercm2s1(long double d) { return d; }

// Luminosity reference
constexpr double operator "" _erg(long double d) { return d; }

#endif
