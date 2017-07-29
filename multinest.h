#ifndef MULTINEST_H
#define MULTINEST_H

//#include <string>
#include <cstring>

namespace nestrun
{
extern "C" {
		void NESTRUN(int &IS, int &mmodal, int &ceff, int &nlive, double &tol, double &efr, int &ndims,
			int &nPar, int &nClsPar, int &maxModes, int &updInt, double &Ztol, char *root, int &seed,
			int *pWrap, int &fb, int &resume, int &outfile, int &initMPI, double &logZero, int &maxiter,
			void (*Loglike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),
			void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void *),
			void *context, int &root_len);
	}

}

class MultiNestFit
{
private:
	static void _loglike(double *Cube, int &n_dim, int &n_par, double &lnew, void *context);
	static void _dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context);
	
protected:
	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 0.8;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int* pWrap;				// which parameters to have periodic boundary conditions?
	
	char root[100]; 		// root for output files
	int root_len = 0;
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 0;				// write output files?
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	const int n_Par;
	virtual double loglike(double *Cube) = 0;
	virtual double dumper(int &nSamples, int &nlive, double** postdist, double** pLivePts, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr) = 0;

public:
	std::string output = "";
	
	MultiNestFit(int n_Par);
	~MultiNestFit();
	MultiNestFit(const MultiNestFit* mnf) = delete;
	MultiNestFit& operator =(const MultiNestFit* mnf) = delete;
	
	void Run();
	
	
};

MultiNestFit::MultiNestFit(int n_Par) : n_Par(n_Par)
{
	pWrap = new int[n_Par];
	for(int i = 0; i < n_Par; i++) pWrap[i] = 0;
}

MultiNestFit::~MultiNestFit() 
{
	delete []pWrap;
}

void MultiNestFit::_loglike(double *Cube, int &n_dim, int &n_par, double &lnew, void *context)
{
	lnew = ((MultiNestFit*)context)->loglike(Cube);
}


void MultiNestFit::_dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// Convert Arrays
	double postdist[nSamples][nPar + 2];
	for(int i = 0; i < nPar + 2; i++ )
		for( int j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( int i = 0; i < nPar + 1; i++ )
		for( int j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
	((MultiNestFit*)context)->dumper(nSamples, nlive,(double**) postdist, (double**) pLivePts, maxLogLike, logZ, INSlogZ, logZerr);
}

void MultiNestFit::Run()
{
	std::strcpy(root, output.c_str());
	root_len = output.length();
	int nPar = n_Par;
	int ndims = n_Par;					// dimensionality (no. of free parameters)
	int nClsPar = n_Par;				// no. of parameters to do mode separation on
	void *context = this;				// not required by MultiNest, any additional information user wants to pass
	// calling MultiNest

	nestrun::NESTRUN(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
					logZero, maxiter, _loglike, _dumper, context, root_len);
}

#endif
