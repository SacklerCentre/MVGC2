#include "LZc_matlab.h"

#include <math.h>
#include <string.h>

size_t maxn(const double* const p, const size_t n)
{
	size_t nmax = 0;
	for (size_t i=0; i<n; ++i) {
		if (p[i] > n) nmax = (double)(p[i]);
	}
	return nmax;
}

// Main function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray* const pistr = prhs[0];                              // pointer to input string or cell string
	const mxArray* const pn = prhs[1];                                 // pointer to input string size(s)
	const size_t nistr = mxGetNumberOfElements(pn);                    // number of strings
	GHashTable* const dict = g_hash_table_new(g_str_hash,g_str_equal); // create dictionary (implemented as a hash set)
	const size_t n = maxn(mxGetPr(pn),nistr);                          // maximum of string lengths
	char* const word = mxMalloc(maxwordlen(n)+1);                      // allocate working word cstring buffer
	for (size_t i=0; i<nistr; ++i) {                                   // for each input string
		char* const istr = mxArrayToString(mxGetCell(pistr,i));        // create input cstring
		LZc(dict,istr,word);                                           // LZ algorithm: build the dictionary
		mxFree(istr);                                                  // deallocate cstring
	}
	mxFree(word);                                                      // deallocate word buffer
	plhs[0] = mxCreateDoubleScalar(g_hash_table_size(dict));           // output LZ complexity (size of dictionary)
	if (nlhs > 1) plhs[1] = hash_table_to_cvec(dict);                  // optionally output dictionary as cell vector of strings
	free_hash_table(dict);                                             // deallocate dictionary
}
