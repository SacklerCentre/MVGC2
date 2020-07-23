#include "mex.h"

/* See tooti.f

mex -R2018a blas_dgemm.c -lblas

*/

typedef long int fint;

void dgemm_(
	char*   TRANSA,
	char*   TRANSB,
	fint*    M,
	fint*    N,
	fint*    K,
	double* ALPHA,
	double* A,
	fint*    LDA,
	double* B,
	fint*    LDB,
	double* BETA,
	double* C,
	fint*    LDC
);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char TRANSA = 'N';
    char TRANSB = 'N';

	double* A   = mxGetDoubles(prhs[0]);
	double* B   = mxGetDoubles(prhs[1]);

	fint M = mxGetM(prhs[0]); /* rows of A */
	fint N = mxGetN(prhs[1]); /* cols of B */
	fint K = mxGetN(prhs[0]); /* cols of A */
	fint L = mxGetM(prhs[1]); /* rows of B */

	double* C = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL));

	fint LDA = M;
	fint LDB = K;
	fint LDC = M;

	double ALPHA = 1.0;

	double BETA = 0.0;

mexPrintf("M = %d\n",M);
mexPrintf("N = %d\n",N);
mexPrintf("K = %d\n",K);
mexPrintf("L = %d\n",L);
mexPrintf("\n");
mexPrintf("LDA = %d\n",LDA);
mexPrintf("LDB = %d\n",LDB);
mexPrintf("LDC = %d\n",LDC);

mexPrintf("*** dgemm in\n",N);
	dgemm_(
		&TRANSA,
		&TRANSB,
		&M,
		&N,
		&K,
		&ALPHA,
		A,
		&LDA,
		B,
		&LDB,
		&BETA,
		C,
		&LDC
	);
mexPrintf("*** dgemm out\n",N);
}
