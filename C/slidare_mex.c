// Use Makefile in this directory to build .mex

#include "mex.h"
#include <string.h> // for memcpy

#define UNUSED __attribute__ ((unused))

typedef long fint; // carefull - assumes C "long" matches SLICOT FORTRAN "INTEGER"!
typedef long flog; // carefull - assumes C "long" matches SLICOT FORTRAN "LOGICAL"!

#define MAX(x1,x2) ((x1) >= (x2) ? (x1) : (x2))
#define MAX4(x1,x2,x3,x4) (MAX(MAX(x1,x2),MAX(x3,x4)))

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
	const fint n2 = 2*n;

	const fint lda = n;
	const fint ldb = n;
	const fint ldq = n;
	const fint ldr = m;
	const fint ldl = n;
	const fint ldx = n;
	const fint lds = n2+m;
	const fint ldt = n2+m;
	const fint ldu = n2;

	// Input matrices

	double* const a = mxGetDoubles(prhs[0]);
	double* const b = mxGetDoubles(prhs[1]);
	double* const q = mxGetDoubles(prhs[2]);
	double* const r = mxGetDoubles(prhs[3]);
	double* const l = mxGetDoubles(prhs[4]);

	// Output matrices

	double* const x      = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)n, mxREAL));
	double* const alphar = mxGetDoubles(plhs[4] = mxCreateDoubleMatrix((mwSize)n2,1,         mxREAL));
	double* const alphai = mxGetDoubles(plhs[5] = mxCreateDoubleMatrix((mwSize)n2,1,         mxREAL));
	double* const beta   = mxGetDoubles(plhs[6] = mxCreateDoubleMatrix((mwSize)n2,1,         mxREAL));

	fint   info;
	double rcond;

	// Although we could, we don't return these; we treat them as workspace matrices

	double* const s = mxCalloc((size_t)(lds*lds),sizeof(double));
	double* const t = mxCalloc((size_t)(ldt*n2), sizeof(double));
	double* const u = mxCalloc((size_t)(ldu*n2), sizeof(double));

	// Workspace matrices

	const fint ldw = (fint)mxGetScalar(prhs[5]);

	const fint liwork = MAX(m,n2);
	const fint ldwork = (ldw > 0 ? ldw : MAX4(7*(n2+1)+16,16*n,n2+m,3*m));
	const fint lbwork = n2;

	fint*   const iwork = mxCalloc((size_t)liwork,sizeof(fint));
	double* const dwork = mxCalloc((size_t)ldwork,sizeof(double));
	flog*   const bwork = mxCalloc((size_t)lbwork,sizeof(flog));

//	mexPrintf("liwork = %ld\n",liwork);
//	mexPrintf("ldwork = %ld\n",ldwork);
//	mexPrintf("lbwork = %ld\n",lbwork);

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

	// Deallocate workspace matrices

	mxFree(bwork);
	mxFree(dwork);
	mxFree(iwork);
	mxFree(u);
	mxFree(t);
	mxFree(s);
}
