#include "mex.h"
#include <string.h> // for memcpy

/*

mex slicot_dlyap.c -lslicot

*/

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const fint    n      = mxGetM(prhs[0]);
	const char    dico   = 'D';
	const char    job    = 'X';
	const char    fact   = 'N';
	const char    trana  = 'N';
	double* const a      = mxGetPr(prhs[0]);
mexPrintf("a =\n\n");
for (fint i=0; i<n; ++i) {
	for (fint j=0; j<n; ++j) mexPrintf("   % 6.4f",a[i+n*j]);
	mexPrintf("\n");
}
mexPrintf("\n");
	const fint    lda    = n;
	double* const u      = mxCalloc((size_t)(n*n),sizeof(double));
	const fint    ldu    = n;
	double* const c      = mxGetPr(prhs[1]);
mexPrintf("c =\n\n");
for (fint i=0; i<n; ++i) {
	for (fint j=0; j<n; ++j) mexPrintf("   % 6.4f",c[i+n*j]);
	mexPrintf("\n");
}
mexPrintf("\n");
	const fint    ldc    = n;
	const double  scale  = 1.0;
	const double  sep    = 0.0;
	const double  ferr   = 0.0;
	double* const wr     = mxCalloc((size_t)n,sizeof(double));
	double* const wi     = mxCalloc((size_t)n,sizeof(double));
	const fint    iwork  = 0;
	const fint    ldwork = 2*n*n+3*n;
	double* const dwork  = mxCalloc((size_t)ldwork,sizeof(double));
	fint          info   = 0;

mexPrintf("*** sb03md in\n");
	sb03md_(
		&dico,
		&job,
		&fact,
		&trana,
		&n,
		a,
		&lda,
		u,
		&ldu,
		c,
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
mexPrintf("*** sb03md out\n");

	mxFree(dwork);
	mxFree(wi);
	mxFree(wr);
	mxFree(u);

	double* const x = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL));
	memcpy(x,c,n*n*sizeof(double));
}
