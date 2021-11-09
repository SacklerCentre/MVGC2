// Use Makefile in this directory to build .mex

#include "mex.h"
#include <string.h> // for memcpy

#define UNUSED __attribute__ ((unused))

typedef long fint;

fint sb03md_(
	const char*   const dico,
	const char*   const job,
	const char*   const fact,
	const char*   const trana,
	const fint*   const n,
	      double* const a,
	const fint*   const lda,
	      double* const u,
	const fint*   const ldu,
	      double* const c,
	const fint*   const ldc,
	const double* const scale,
	const double* const sep,
	const double* const ferr,
	      double* const wr,
	      double* const wi,
	const fint*   const iwork,
	      double* const dwork,
	const fint*   const ldwork,
	      fint*   const info
);

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	const mwSize  N      = mxGetM(prhs[0]);
	const mwSize  NSQ    = N*N;
	const fint    n      = (fint)N;
	const char    dico   = 'D';
	const char    job    = 'X';
	const char    fact   = 'N';
	const char    trana  = 'T';

	const double* const a = mxGetDoubles(prhs[0]);
	const fint    lda     = n;
	double* const u       = mxCalloc(NSQ,sizeof(double));
	const fint    ldu     = n;
	const double* const c = mxGetDoubles(prhs[1]);
	const fint    ldc     = n;
	const double  scale   = 1.0;
	const double  sep     = 0.0;
	const double  ferr    = 0.0;
	double* const wr      = mxCalloc(N,sizeof(double));
	double* const wi      = mxCalloc(N,sizeof(double));

	const fint    iwork  = 0;
	const fint    ldwork = 2*(fint)(NSQ+3*N); // reasonable?
	double* const dwork  = mxCalloc((mwSize)ldwork,sizeof(double));
	fint          info   = 0;

	// Note that SLICOT SB03MD overwrites A, and also C (the result is returned in C)

	// Take a copy of a
	double* const A = mxCalloc(NSQ,sizeof(double));
	memcpy(A,a,NSQ*sizeof(double));

	// Take a copy of c, and assign as first output
	double* const C = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL));
	memcpy(C,c,NSQ*sizeof(double));

	// mexPrintf("*** sb03md in\n");
	sb03md_(
		&dico,
		&job,
		&fact,
		&trana,
		&n,
		A,
		&lda,
		u,
		&ldu,
		C,
		&ldc,
		&scale,
		&sep,
		&ferr,
		wr,
		wi,
		&iwork,
		dwork,
		&ldwork,
		&info
	);
	// mexPrintf("*** sb03md out\n");

	mxFree(A);
	mxFree(dwork);
	mxFree(wi);
	mxFree(wr);
	mxFree(u);
}
