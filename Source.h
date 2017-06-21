#ifndef SOURCE_H
#define SOURCE_H

#include <iostream>
#include <functional>
#include <tuple>
#include <utility>
#include <memory>
#include "Constants.h"
#include "CosmologyModel.h"


typedef std::pair<double, double> Bounds ;			// To store Integration Bounds (left, right)
typedef std::tuple<double, double, double> Bin;  	// To store energy bins  (left, mid, right)

class DGRBSource
{
	friend class Benchmark;
public:
	std::string Name;
	DGRBSource(std::shared_ptr<CosmologyModel> _CM) : DGRBSource(_CM, std::string("")) {} 
	DGRBSource(std::shared_ptr<CosmologyModel> _CM, std::string _name) : Name(_name), CM(_CM) {}
	
	std::shared_ptr<CosmologyModel> CM;
	
	std::function<double(const double *)> Intensity;
	std::function<double(const double *)> Autocorrelation;
	
	virtual ~DGRBSource() {}
};



#endif
