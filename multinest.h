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

//#endif // ifdef __cplusplus

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include "Constants.h"


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

struct MultiNestFitData
{
	//MultiNestFitDump(int nSamples, int nLivePoints, int nPar, double
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

class MultiNestFit
{
private:
	static void _loglike(double *Cube, int *n_dim, int *n_par, double *lnew, void *context);
	static void _dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, void *context);
	
protected:
	unsigned int nParameters;
	
	virtual double loglike(double *Cube, int *npar) = 0;
	virtual void dumper(MultiNestFitData& data) = 0;

public:
	MultiNestFitOptions options;
	MultiNestFitData data;
	
	MultiNestFit(std::string output);
	//~MultiNestFit();
	//MultiNestFit(const MultiNestFit* mnf) = delete;
	//MultiNestFit& operator =(const MultiNestFit* mnf) = delete;
	
	MultiNestFitData& Run();
	
	
};

MultiNestFit::MultiNestFit( std::string output = "") 
{	
	options.IS = 1;					// do Nested Importance Sampling?
	options.mmodal = 0;					// do mode separation?
	options.ceff = 0;					// run in constant efficiency mode?
	options.nlive = 1000;				// number of live points
	options.efr = 0.8;				// set the required efficiency
	options.tol = 0.5;				// tol, defines the stopping criteria
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


void MultiNestFit::_loglike(double *Cube, int *n_dim, int *n_par, double *lnew, void *context)
{
	//for(unsigned int i = 0; i < ((MultiNestFit*)context)->parameters.size(); i++) ((MultiNestFit*)context)->parameters[i].lastValue = Cube[i];
	*lnew = ((MultiNestFit*)context)->loglike(Cube, n_par);
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

MultiNestFitData& MultiNestFit::Run()
{
	char root[100];
	std::strcpy(root, options.root.c_str());
	for(int i = std::strlen(root); i < 100; i++) root[i] = ' ';
	
	int nPar = (int)nParameters;
	int ndims = nParameters;					// dimensionality (no. of free parameters)
	int nClsPar = nParameters;				// no. of parameters to do mode separation on
	int pWrap[nParameters];
	for(unsigned int i = 0; i < nParameters; i++) pWrap[0]=0; //pWrap[i] = (int)parameters[i].periodicBoundaryConditions;
	void *context = this;				// not required by MultiNest, any additional information user wants to pass
	
	data.lastParameters.resize(nPar);
	// calling MultiNest

	//nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	//				logZero, maxiter, _loglike, _dumper, context);
	for (int i = strlen(root); i < 100; i++) root[i] = ' ';

     NESTRUN(&options.IS, &options.mmodal, &options.ceff, &options.nlive, &options.tol, &options.efr, &ndims, &nPar, &nClsPar, &options.maxModes, &options.updInt, &options.Ztol,
        root, &options.seed, pWrap, &options.fb, &options.resume, &options.outfile, &options.initMPI, &options.logZero, &options.maxiter, &_loglike, &_dumper, context);
     
    return data;
}

#endif
