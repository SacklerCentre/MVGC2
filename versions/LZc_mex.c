#include "mex.h"
#include <math.h>
#include <string.h>
#include <gmodule.h>

/*
	Calculate LZ complexity, and (optionally) return the dictionary. Call in Matlab as (see LZc.m)

		[n,dict] = LZc(s,true);

	NOTE: This is C99, and needs Glib "gmodule" for GHashTable.

	To compile on Linux with GCC (static linkage for portability):

		mex -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Wall -Werror -O3 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -static" -lglib-2.0 C/LZc_mex.c

	This should also work on Windows with MinGW-w64 (https://mingw-w64.org/doku.php),
	(note: MS Visual Studio doesn't do C99), and probably fine on Mac, but not tested.

		*** WARNING: out-of-memory not checked (TODO)! ***

	For 1,000,000 character random binary input strings, this function is ~100 x faster than
	an equivalent Matlab implementation using containers.Map (see LZc.m). E.g.,

		Matlab  (using containers.Map) : 7.11 seconds.
		LZc_mex (using GHashTable)     : 0.07 seconds

*/

void free_hash_table(GHashTable* dict)
{
	GHashTableIter iter;
	gpointer key;
	g_hash_table_iter_init(&iter,dict);
	while (g_hash_table_iter_next (&iter,&key,NULL)) free(key);
	g_hash_table_destroy(dict);
}

mxArray* const hash_table_to_cvec(GHashTable* dict)
{
	mxArray* const cvec = mxCreateCellMatrix(g_hash_table_size(dict),1);
	size_t i = 0;
	GHashTableIter iter;
	gpointer key;
	g_hash_table_iter_init(&iter,dict);
	while (g_hash_table_iter_next (&iter,&key,NULL)) mxSetCell(cvec,i++,mxCreateString(key));
	return cvec;
}

const size_t maxwordlen(const size_t n)
{
	return floor((sqrt(8*n+1)-1.0)/2.0); // worst-case scenario
}

void LZc(GHashTable* const dict, const char* const istr, char* const word)
{
	// word MUST be large enough! See maxwordlen()
	char* w;                                           // pointer to end of current word (null terminator)
	*(w = word) = '\0';                                // initialise to empty word
	for (const char* c = istr; *c; ++c) {              // traverse input string (terminating condition equiv to *c == '\0' !!!)
		*w = *c; *(++w) = '\0';                        // append next input character to word
		if (g_hash_table_lookup(dict,word) == NULL) {  // if word not in dictionary,
			g_hash_table_add(dict,strdup(word));       // add word to dictionary (treated as a hash set!)
			*(w = word) = '\0';                        // and reinitialise to empty word
		}
	}
}

const size_t maxn(const double* const p, const size_t n)
{
	size_t mxn = 0;
	for (size_t i=0; i<n; ++i) {
		if (p[i] > n) mxn = p[i];
	}
	return mxn;
}

// Main function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray* const pistr = prhs[0];                                 // pointer to input string or cell string
	const mxArray* const pn = prhs[1];                                    // pointer to input string size(s)
	const size_t nistr = mxGetNumberOfElements(pn);                       // number of strings
	GHashTable* const dict = g_hash_table_new(g_str_hash,g_str_equal);    // create dictionary (implemented as a hash set)
	if (nistr > 1) {                                                      // if multiple strings (cell string)
		const size_t n = maxn(mxGetPr(pn),nistr);                         // maximum of string lengths
		char* const word = mxMalloc(maxwordlen(n)+1);                     // allocate working word cstring buffer
		for (size_t i=0; i<nistr; ++i) {                                  // for each input string
			char* const istr = mxArrayToString(mxGetCell(pistr,i));       // create input cstring
			LZc(dict,istr,word);                                          // LZ algorithm: build the dictionary
			mxFree(istr);                                                 // deallocate cstring
		}
		mxFree(word);                                                     // deallocate word buffer
	}
	else {                                                                // else just a string
		const size_t n = mxGetScalar(pn);                                 // string length
		char* const word = mxMalloc(maxwordlen(n)+1);                     // allocate working word cstring buffer
		char* const istr = mxArrayToString(pistr);                        // create input cstring
		LZc(dict,istr,word);                                              // LZ algorithm: build the dictionary
		mxFree(istr);                                                     // deallocate cstring
		mxFree(word);                                                     // deallocate word buffer
	}
	plhs[0] = mxCreateDoubleScalar(g_hash_table_size(dict));              // output LZ complexity (size of dictionary)
	if (nlhs > 1) plhs[1] = hash_table_to_cvec(dict);                     // optionally output dictionary as cell vector of strings
	free_hash_table(dict);                                                // deallocate dictionary
}
