// Use Makefile in this directory to build .mex

#include "mex.h"
#include <string.h> // for memcpy

#define UNUSED __attribute__ ((unused))

typedef long fint; // carefull - assumes C "long" matches SLICOT FORTRAN "INTEGER"!
typedef long flog; // carefull - assumes C "long" matches SLICOT FORTRAN "LOGICAL"!

fint sb02od_(
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

	const fint n  = (fint)mxGetM(prhs[0]);
	const fint m  = (fint)mxGetM(prhs[3]);
	const fint p  = 0; // not used for fact == 'N'

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

	double* const a = mxGetDoubles(prhs[0]);
	double* const b = mxGetDoubles(prhs[1]);
	double* const q = mxGetDoubles(prhs[2]);
	double* const r = mxGetDoubles(prhs[3]);
	double* const l = mxGetDoubles(prhs[4]);

	// Output matrices

	double* const x = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix((mwSize)n,(mwSize)n,mxREAL));

	fint   info;
	double rcond;

	// Actually we don't return these; we treat them as workspace matrices, preallocated

	double* const alphar = mxGetDoubles(prhs[5]);
	double* const alphai = mxGetDoubles(prhs[6]);
	double* const beta   = mxGetDoubles(prhs[7]);
	double* const s      = mxGetDoubles(prhs[8]);
	double* const t      = mxGetDoubles(prhs[9]);
	double* const u      = mxGetDoubles(prhs[10]);

	// Workspace matrices, also preallocated

	const fint ldwork = (fint)mxGetNumberOfElements(prhs[12]);

//	mexPrintf("ldwork = %ld\n",ldwork);

	fint*   const iwork = mxGetInt64s (prhs[11]); // carefull - assumes C "long" is Matlab "int64"!
	double* const dwork = mxGetDoubles(prhs[12]);
	flog*   const bwork = mxGetInt64s (prhs[13]); // carefull - assumes C "long" is Matlab "int64"!

	// Other parameters

	const double tol = 0.0; // default tolerance

//	mexPrintf("*** sb02od in\n");
	sb02od_(
		&dico,
		&jobb,
		&fact,
		&uplo,
		&jobl,
		&sort,
		&n,
		&m,
		&p,
	    a,
		&lda,
	    b,
		&ldb,
		q,
		&ldq,
		r,
		&ldr,
		l,
		&ldl,
		&rcond,
		x,
		&ldx,
		alphar,
		alphai,
		beta,
		s,
		&lds,
		t,
		&ldt,
		u,
		&ldu,
		&tol,
		iwork,
		dwork,
		&ldwork,
		bwork,
		&info
	);
//	mexPrintf("*** sb02od out\n");

	plhs[1] = mxCreateDoubleScalar((double)info);
	plhs[2] = mxCreateDoubleScalar(rcond);
	plhs[3] = mxCreateDoubleScalar(dwork[0]); // if info == 0, this is the optimal ldwork
}
