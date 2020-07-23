#include "mex.h"
#include <string.h> // for memcpy

/*

mex slicot_dlyap.c -lslicot

*/

#define UNUSED __attribute__ ((unused))

typedef mwSignedIndex fint;
typedef ??? flog;

/*
  SB02OD

  To solve for X either the continuous-time algebraic Riccati
  equation
                           -1
     Q + A'X + XA - (L+XB)R  (L+XB)' = 0                       (1)

  or the discrete-time algebraic Riccati equation
                                  -1
     X = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q              (2)

  where A, B, Q, R, and L are N-by-N, N-by-M, N-by-N, M-by-M and
  N-by-M matrices, respectively, such that Q = C'C, R = D'D and
  L = C'D; X is an N-by-N symmetric matrix.
  The routine also returns the computed values of the closed-loop
  spectrum of the system, i.e., the stable eigenvalues lambda(1),
  ..., lambda(N) of the corresponding Hamiltonian or symplectic
  pencil, in the continuous-time case or discrete-time case,
  respectively.
                           -1
  Optionally, matrix G = BR  B' may be given instead of B and R.
  Other options include the case with Q and/or R given in a
  factored form, Q = C'C, R = D'D, and with L a zero matrix.

  The routine uses the method of deflating subspaces, based on
  reordering the eigenvalues in a generalized Schur matrix pair.
  A standard eigenproblem is solved in the continuous-time case
  if G is given.
*/

fint sb02od_(
};

fint sb03md_(
	const char*   const dico,
	const char*   const jobb,
	const char*   const fact,
	const char*   const uplo,
	const char*   const jobl,
	const char*   const sort,
	const fint*   const n,
	const fint*   const m,
	const fint*   const p,
	      double* const a,
	const fint*   const lda,
	      double* const b,
	const fint*   const ldb,
	      double* const q,
	const fint*   const ldq,
	      double* const r,
	const fint*   const ldr,
	      double* const l,
	const fint*   const ldl,
	const double* const rcond,
	      double* const x,
	const fint*   const ldx,
	      double* const alphar,
	      double* const alphai,
	      double* const beta,
	      double* const s,
	const fint*   const lds,
	      double* const t,
	const fint*   const ldt,
	      double* const u,
	const fint*   const ldu,
	const double* const tol,
	      fint*   const iwork,
	      double* const dwork,
	const fint*   const ldwork,
	      flog*   const bwork,
	      fint*   const info
);

// Parameters: A, B, Q, R, and L are N-by-N, N-by-M, N-by-N, M-by-M and N-by-M matrices, respectively

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// Job constants

	const char dico = 'D';
	const char jobb = 'B';
	const char fact = 'N';
	const char uplo = 'L';
	const char jobl = 'N';
	const char sort = 'S';

	// Array sizes

	const fint n = (fint)mxGetM(prhs[0]);
	const fint m = (fint)mxGetM(prhs[3]);
	const fint p = 0; // not used for fact == 'N'

	const size_t N  = (size_t)n;
	const size_t N2 = 2*N;
	const size_t M  = (size_t)m;
	const size_t NM2 = N2+M;

	const fint lda = n;
	const fint ldb = n;
	const fint ldq = n;
	const fint ldr = m;
	const fint ldl = n;
	const fint ldx = n;
	const fint lds = 2*n+m;
	const fint ldt = 2*n+m;
	const fint ldu = 2*n;

	// Input matrices

	double* const a = mxGetPr(prhs[0]);
	double* const b = mxGetPr(prhs[1]);
	double* const q = mxGetPr(prhs[2]);
	double* const r = mxGetPr(prhs[3]);
	double* const l = mxGetPr(prhs[4]);

	// Output matrices

	double* const x      = mxGetPr(plhs[0] = mxCreateDoubleMatrix(N, N,   mxREAL));
	double* const rcond  = mxGetPr(plhs[1] = mxCreateDoubleMatrix(1, 1,   mxREAL));
	double* const alphar = mxGetPr(plhs[2] = mxCreateDoubleMatrix(N2,N2,  mxREAL));
	double* const alphai = mxGetPr(plhs[3] = mxCreateDoubleMatrix(N2,N2,  mxREAL));
	double* const beta   = mxGetPr(plhs[4] = mxCreateDoubleMatrix(N2,N2,  mxREAL));
	double* const s      = mxGetPr(plhs[5] = mxCreateDoubleMatrix(N2,NM2, mxREAL));
	double* const t      = mxGetPr(plhs[6] = mxCreateDoubleMatrix(N2,NM2, mxREAL));
	double* const u      = mxGetPr(plhs[6] = mxCreateDoubleMatrix(N2,N2,  mxREAL));

	// Workspace matrices

	fint*   const iwork = mxCalloc(MAX(M,N2),sizeof(fint));
	double* const dwork = mxCalloc(MAX(7*(2*N+1)+16,16*N,2*N+M,3*M),sizeof(double));
	double* const bwork = mxCalloc(N2,sizeof(double));


	const   fint   n      = (fint)mxGetM(prhs[0]);
	const   size_t N      = (size_t)n;
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
