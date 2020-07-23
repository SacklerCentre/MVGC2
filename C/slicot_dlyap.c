#include "mex.h"
#include <string.h> // for memcpy

/*

mex slicot_dlyap.c -lslicot

*/

#define UNUSED __attribute__ ((unused))

typedef mwSignedIndex fint;

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
	const   fint   n      = (fint)mxGetM(prhs[0]);
	const   size_t N      = (size_t)n;
	const   char   dico   = 'D';
	const   char   job    = 'X';
	const   char   fact   = 'N';
	const   char   trana  = 'T';
	double* const  a      = mxGetPr(prhs[0]);
	const   fint   lda    = n;
	double* const  u      = mxCalloc((size_t)(n*n),sizeof(double));
	const   fint   ldu    = n;
	double* const  c      = mxGetPr(prhs[1]);
	const   fint   ldc    = n;
	const   double scale  = 1.0;
	const   double sep    = 0.0;
	const   double ferr   = 0.0;
	double* const  wr     = mxCalloc(N,sizeof(double));
	double* const  wi     = mxCalloc(N,sizeof(double));
	const   fint   iwork  = 0;
	const   fint   ldwork = 2*n*n+3*n;
	double* const  dwork  = mxCalloc((size_t)ldwork,sizeof(double));
	fint           info   = 0;
	double* const  A      = mxCalloc(N*N,sizeof(double));
	double* const  C      = mxGetPr(plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL));

	memcpy(A,a,N*N*sizeof(double));
	memcpy(C,c,N*N*sizeof(double));

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

	mxFree(dwork);
	mxFree(wi);
	mxFree(wr);
	mxFree(u);
	mxFree(A);
}
