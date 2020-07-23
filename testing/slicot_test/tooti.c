#include "mex.h"

/* See tooti.f

mex LDFLAGS="\$LDFLAGS -L." mxtooti.c -ltooti

*/

void tooti_(double* const,const double* const);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const double* const x = mxGetPr(prhs[0]);

	double* const y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL));

	tooti_(y,x);
}
