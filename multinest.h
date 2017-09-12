#ifndef MULTINEST_H
#define MULTINEST_H

#ifdef __INTEL_COMPILER 			// if the MultiNest library was compiled with ifort
//       #define NESTRUN nested_mp_nestrun_
       #define NESTRUN __nested_MOD_nestrun         
#elif defined __GNUC__ 				// if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) 
#endif



/***************************************** C Interface to MultiNest **************************************************/
extern "C" {
extern void NESTRUN(int *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, 
char *, int *, int *, int *, int *, int *, int *, double *, int *, void (*Loglike)(double *, int *, int *, 
double *, void *), void (*dumper)(int *, int *, int *, double **, double **, double **, double *, 
double *, double *, double *, void *), void *context);
}


/***********************************************************************************************************************/

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <memory>
#include <cassert>
#include "Constants.h"

// A struct to hold the options, so the class namespace isn't filled with crap
struct MultiNestFitOptions
{
	int IS;  		// Nested Importance Sampling
	int mmodal; 	// mode separation
	int ceff; 		// constant efficiency mode?
	int nlive;		// number of live points
	double efr;		// required efficiency
	double tol; 	// stopping criteria
	int updInt;		// posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol;	// modes with logZ , Ztol are ignored
	int maxModes;	//expected max number of nodes
	
	std::string root; 	// dir for output
	int seed;		// seed for random numbers
	int fb;			// feedback on std
	int resume;		// resume from a previous job
	int outfile;	// write output
	int initMPI;	// init MPI routines
	double logZero;	// points with loglike < logZero will be ignored
	int maxiter;	// max no. of iteration / non-positive means infinte
};

// Collects data for the dumper 
// Should be changed to hold all past live points
// But wasn't needed so far
struct MultiNestFitData
{
	//double** posteriorDistribution;
	std::vector<double> lastParameters;
	double logLike;
	double maxLogLike;
	double logZ;
	double INSlogZ;
	double logZerr;
	int nLivePoints;
	int nSamples;
	 int nPar;
};

// A small parameter class, which takes the internal multinest parameter and scales it linearly inside its bounds
class MultiNestParameter final	// is final atm, so that the std::vector<MultiNestParameter> in the MultiNestFit class does not fall victim to object splicing 
{
friend class MultiNestFit;
protected:
	double internalMultinestParameter;	// between 0 and 1
	double getValue(double internalMultinestParameter)
	{
		this->internalMultinestParameter = internalMultinestParameter;
		return getValue();
	}
public:
	std::string name;
	Bounds bounds;
	bool periodicBoundaryCondition;
	double getValue()
	{
		return bounds.first + (bounds.second - bounds.first)*internalMultinestParameter;
	}
	double operator()()	{return getValue(); }	// syntactic sugar
	// standard constructors will work
	
	MultiNestParameter(Bounds bounds, std::string name = "", bool periodicBoundaryCondition = false) : name(name), bounds(bounds), periodicBoundaryCondition(periodicBoundaryCondition)
	{
		internalMultinestParameter = 0;
	}
	void print()
	{
		std::cout << name << ": " << getValue() << std::endl;
	}
};


// This class uses the *context pointer in multinest to call its own loglike and dumper function
// So that a fit class can hold data and use it easily without using too much pointer magic
// Other classes should inherit from this
class MultiNestFit
{
private:
	static void _loglike(double *Cube, int *n_dim, int *n_par, double *lnew, void *context);
	static void _dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, void *context);
	
protected:
	std::vector<MultiNestParameter> Parameters;		// The list of parameters can be used, but doesn't have to be
													// The size of this determines nPar passed to multinest though
													// And on _loglike call the list is updated
													// So if not used, the Run, _loglike, and _dumper functions should be overwritten
	virtual double loglike() = 0;	// child classes need to implement these
	virtual void dumper(MultiNestFitData& data) = 0;

public:
	MultiNestFitOptions options;
	MultiNestFitData data;
	
	MultiNestFit(std::string output);
	MultiNestFit(const MultiNestFitOptions& mnfo);
	//~MultiNestFit();		// The default constructors should work
	//MultiNestFit(const MultiNestFit* mnf) = delete;
	//MultiNestFit& operator =(const MultiNestFit* mnf) = delete;
	
	MultiNestFitData& Run();
	
	
};

MultiNestFit::MultiNestFit( std::string output = "") // setting some standards
{	
	options.IS = 1;					// do Nested Importance Sampling?
	options.mmodal = 0;					// do mode separation?
	options.ceff = 0;					// run in constant efficiency mode?
	options.nlive = 1000;				// number of live points
	options.efr = 0.3;				// set the required efficiency
	options.tol = 0.1;				// tol, defines the stopping criteria
	options.updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	options.Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	options.maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	
	options.seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	options.fb = 1;					// need feedback on standard output?
	options.resume = 0;					// resume from a previous job?
	options.outfile = 0;				// write output files?
	options.initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	options.logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	options.maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
							
	if(!output.empty())
	{
		this->options.root = output;
		this->options.outfile = 1;
	}
}

MultiNestFit::MultiNestFit(const MultiNestFitOptions& mnfo) : options(mnfo)
{	}

void MultiNestFit::_loglike(double *Cube, int *n_dim, int *n_par, double *lnew, void *context)	// uses the *context pointer to call the MultiNestFit::loglike function
{
	auto mnf = ((MultiNestFit*)context);
	assert((int)mnf->Parameters.size() == *n_dim);	// this needs to be changed if the parameter lists isn't used
	for(unsigned int i = 0; i < mnf->Parameters.size(); i++)  Cube[i] = mnf->Parameters[i].getValue(Cube[i]);	// save parameters in class and return the rescaled version
	*lnew = mnf->loglike();	// actually calculate the loglike
}


void MultiNestFit::_dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, void *context)
{	
	MultiNestFit* mnf = ((MultiNestFit*)context);
	for(int i = 0; i < *nPar; i++) mnf->data.lastParameters[i] = physLive[0][i * (*nlive)];
	mnf->data.logLike = physLive[0][*nPar * (*nlive)];
	
	mnf->data.maxLogLike = *maxLogLike;
	mnf->data.logZ = *logZ;
	mnf->data.INSlogZ = *INSlogZ;
	mnf->data.logZerr = *logZerr;
	mnf->data.nLivePoints = *nlive;
	mnf->data.nSamples = *nSamples;
	mnf->data.nPar = *nPar;
	
	mnf->dumper(mnf->data);
}

// This starts the multinest fitting with the specified parameters
MultiNestFitData& MultiNestFit::Run()
{
	char root[100];
	std::strcpy(root, options.root.c_str());
	for(int i = std::strlen(root); i < 100; i++) root[i] = ' ';		// prepare output string as multinest wants it
	
	int nPar = Parameters.size();
	int ndims = nPar;					// dimensionality (no. of free parameters)
	int nClsPar = nPar;				// no. of parameters to do mode separation on
	int pWrap[nPar];
	for(unsigned int i = 0; i < Parameters.size(); i++) pWrap[i]=(int)Parameters[i].periodicBoundaryCondition; 
	void *context = this;				// This is why we can use the class here
	data.lastParameters.resize(nPar);
	
	// calling MultiNest

     NESTRUN(&options.IS, &options.mmodal, &options.ceff, &options.nlive, &options.tol, &options.efr, &ndims, &nPar, &nClsPar, &options.maxModes, &options.updInt, &options.Ztol,
        root, &options.seed, pWrap, &options.fb, &options.resume, &options.outfile, &options.initMPI, &options.logZero, &options.maxiter, &_loglike, &_dumper, context);
     
    return data;
}

#endif
