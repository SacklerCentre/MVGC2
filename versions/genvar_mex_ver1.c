#include <string.h>
#include "mex.h"

/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling genvar_mex itself! It is designed to be called via the
 * genvar.m Matlab wrapper function (utils subfolder); all error checking - not
 * to mention documentation - happens there. To compile, see mvgc_makemex.m in
 * the utils subfolder */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double* a = mxGetPr(prhs[0]);
	double* x = mxGetPr(prhs[1]);

	const size_t n  = mxGetM(prhs[1]);
	const size_t m  = mxGetN(prhs[1]);
	const mwSize d  = mxGetNumberOfDimensions(prhs[0]);
	const mxLogical aisvec = *mxGetPr(prhs[2]);
	const mwSize p  = (aisvec ? mxGetNumberOfElements(prhs[0]) : d < 3 ? 1 : mxGetDimensions(prhs[0])[2]);

	double* y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); /* allocate output array y  */

	const size_t nsq = n*n;
    size_t t,k,i,j,toff,koff,lagoff,ioff;
    double yit;

if (aisvec) {
	mexPrintf("true : p = %d\n",p);
}
else {
	mexPrintf("false : p = %d\n",p);
}

	memcpy(y,x,n*m*sizeof(double)); /* copy input x to output y */
	for (t=0; t<p; ++t) {
		toff = n*t;
		for (k=0; k<t; ++k) {
			koff = nsq*k;
			lagoff = toff-n*(k+1); /* lag is k+1 !!! */
			for (i=0; i<n; ++i) {
				yit = 0.0;
				ioff = i+koff;
				for (j=0; j<n; ++j) {
					yit += a[ioff+n*j]*y[j+lagoff];
				}
				y[i+toff] += yit;
			}
		}
	}
	for (; t<m; ++t) {
		toff = n*t;
		for (k=0; k<p; ++k) {
			koff = nsq*k;
			lagoff = toff-n*(k+1); /* lag is k+1 !!! */
			for (i=0; i<n; ++i) {
				yit = 0.0;
				ioff = i+koff;
				for (j=0; j<n; ++j) {
					yit += a[ioff+n*j]*y[j+lagoff];
				}
				y[i+toff] += yit;
			}
		}
	}
}
