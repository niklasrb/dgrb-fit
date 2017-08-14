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


//#include <string>
#include <cstring>

#ifdef __cplusplus

/***************************************** C++ Interface to MultiNest **************************************************/

#include <cstring>
#include <string>

namespace nested
{

	// map the Fortran 90 entry points of libnest3.a to C++ functions

	// module nested, function nestRun maps to nested::run

	// the pass-by-reference nature of most of the Fortran is translated away
	// *apart* from the callbacks. The provided call back functions must still accept 
	// references rather than values. There is also some confusion as to the type
	// of the first argument of LogLike. 
	// Should it be a double * or an farray<double, 1> *? The former seems to 
	// work and is simpler.

	// This structure is reverse engineered from looking 
	// at gfortran stack traces. It is probably wrong
	
	template<typename type, int ndims> class farray_traits;
	
	template<> class farray_traits<double, 1> { public: static const int id = 537; };
	template<> class farray_traits<double, 2> { public: static const int id = 538; };
	template<> class farray_traits<int, 1> { public: static const int id = 265; };
	template<> class farray_traits<int, 2> { public: static const int id = 266; };

	// the extra data for f90 that defines how arrays are arranged.
	template<typename T, int ndim> class farray
	{
		public:
			farray(T *_data, int w, int h = 0) : data(_data), offset(0), type(farray_traits<T, ndim>::id), 
			x_stride(1), x_lbound(1), x_ubound(w), y_stride(w), y_lbound(1), y_ubound(h) {};
			
			T *data;
			int offset;
			int type;
			int x_stride, x_lbound, x_ubound;
			int y_stride, y_lbound, y_ubound;
	};
	
	extern "C" {
		void NESTRUN(int &IS, int &mmodal, int &ceff, int &nlive, double &tol, double &efr, int &ndims,
			int &nPar, int &nClsPar, int &maxModes, int &updInt, double &Ztol, char *root, int &seed,
			int *pWrap, int &fb, int &resume, int &outfile, int &initMPI, double &logZero, int &maxiter,
			void (*Loglike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),
			void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void *),
			void *context, int &root_len);
	}

	static void run(bool IS, bool mmodal, bool ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, int maxModes,
		int updInt, double Ztol, const std::string & root, int seed, int *pWrap, bool fb, bool resume, bool outfile, 
		bool initMPI, double logZero, int maxiter, void (*LogLike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),
		void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void *), void *context)
	{
		char t_root[100];
		std::fill(t_root, t_root + 100, ' ');
		snprintf(t_root, 99, "%s", root.c_str());
		int root_len = strlen(t_root);
		t_root[strlen(t_root)] = ' ';
	
		int t_fb = fb;
		int t_resume = resume;
		int t_outfile = outfile;
		int t_initMPI = initMPI;
		int t_mmodal = mmodal;
		int t_IS = IS;
		int t_ceff = ceff;
		
		NESTRUN(t_IS, t_mmodal, t_ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, t_root, seed, pWrap, t_fb, 
		t_resume, t_outfile, t_initMPI, logZero, maxiter, LogLike, dumper, context, root_len);
	}	
}

/***********************************************************************************************************************/

#else // ifdef __cplusplus

/***************************************** C Interface to MultiNest **************************************************/

extern void NESTRUN(int *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, 
char *, int *, int *, int *, int *, int *, int *, double *, int *, void (*Loglike)(double *, int *, int *, 
double *, void *), void (*dumper)(int *, int *, int *, double **, double **, double **, double *, 
double *, double *, double *, void *), void *context);

void run(int IS, int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, 
int maxModes, int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, 
int initMPI, double logZero, int maxiter, void (*LogLike)(double *, int *, int *, double *, void *), 
void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, double *, void *), 
void *context)
{
	int i;
	for (i = strlen(root); i < 100; i++) root[i] = ' ';

        NESTRUN(&IS, &mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, &maxiter, LogLike, dumper, context);
}

/***********************************************************************************************************************/

#endif // ifdef __cplusplus


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
	
	std::string root; 		// root for output files
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
	//std::string output = "";
	
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
	//std::strcpy(root, output.c_str());
	//root_len = output.length();
	int nPar = n_Par;
	int ndims = n_Par;					// dimensionality (no. of free parameters)
	int nClsPar = n_Par;				// no. of parameters to do mode separation on
	void *context = this;				// not required by MultiNest, any additional information user wants to pass
	// calling MultiNest

	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
					logZero, maxiter, _loglike, _dumper, context);
}

#endif
