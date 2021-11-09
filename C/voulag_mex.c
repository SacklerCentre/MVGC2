// Use Makefile in this directory to build .mex

#include <string.h> // for memcpy
#include "mex.h"
#include "matrix.h"
#include <inttypes.h>

/**************************************************************************
 * WARNING: NO ERROR CHECKING WHATSOEVER - CALL IT RIGHT!!! (see oulag.m) *
 **************************************************************************
 *
 *      X = voulag_mex(Z,d,a,uint64(to-1),uint64(from-1),uint64(lagm));
 *
 **************************************************************************/

// #define ct_assert(e) extern char (*ct_assert(void)) [sizeof(char[1 - 2*!(e)])] // weirs but effective compile-time assert

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	if (sizeof(size_t) != sizeof(uint64_t)) { // ensure Matlab integer type matches size_t !!
		mexErrMsgTxt("Bad array index size; please amend this routine and \"vougc.m\" for your platform!");
	}

	// get inputs
	const double*  const z    = mxGetPr(prhs[0]);
	const double*  const d    = mxGetPr(prhs[1]);
	const double*  const a    = mxGetPr(prhs[2]);
	const size_t*  const to   = (size_t*)mxGetData(prhs[3]);
	const size_t*  const from = (size_t*)mxGetData(prhs[4]);
	const size_t*  const lagm = (size_t*)mxGetData(prhs[5]);

	const size_t n = mxGetM(prhs[0]);                // num variables
	const size_t m = mxGetN(prhs[0]);                // num observations
	const size_t c = mxGetNumberOfElements(prhs[2]); // num connections

	// allocate output variable x and copy it from z
	double* const x = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL));
	memcpy(x,z,n*m*sizeof(double));

    for (size_t t=1;t<m;++t) {
		double* const xt = x+n*t;
		for (size_t i=0;i<n;++i) xt[i] += d[i]*xt[i-n]; // decay loop
		for (size_t k=0;k<c;++k) { // causal loop
			if (lagm[k] <= t) xt[to[k]] += a[k]*xt[from[k]-n*lagm[k]];
		}
	}
}
