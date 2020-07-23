#include <string.h> // for memcpy
#include "mex.h"

/**************************************************************************
 * WARNING: NO ERROR CHECKING WHATSOEVER - CALL IT RIGHT!!! (see oulag.m) *
 **************************************************************************
 *
 *     X = oulag_c(Z,b,A,r);
 *
 */

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	const double* const z = mxGetPr(prhs[0]);
	const double* const b = mxGetPr(prhs[1]);
	const double* const a = mxGetPr(prhs[2]);

	const size_t n = mxGetM(prhs[0]);              // num variables
	const size_t m = mxGetN(prhs[0]);              // num observations
	const size_t r = (size_t)mxGetScalar(prhs[3]); // num lags

	double* const x = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); // allocate output variable x
	memcpy(x,z,n*m*sizeof(double));                                        // and copy it from z

    for (size_t k=r;k<m;++k) {
		for (size_t i=0;i<n;++i) {
            double axi = 0.0;
            for (size_t j=0;j<n;++j) axi += a[i+n*j]*x[j+n*(k-r)];
            x[i+n*k] += b[i]*x[i+n*(k-1)]+axi;
        }
    }
}
