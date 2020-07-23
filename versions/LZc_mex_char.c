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

void LZc(GHashTable* const dict, const char* const istr, char* const word)
{
	char* w;                                           // current word terminator
	*(w = word) = '\0';                                // initialise to empty word
	for (const char* c = istr;*c;++c) {                // traverse input string (terminating condition equiv to *c == '\0' !!!)
		*w = *c; *(++w) = '\0';                        // append next input character to word
		if (g_hash_table_lookup(dict,word) == NULL) {  // if word not in dictionary,
			g_hash_table_add(dict,strdup(word));       // add to dictionary
			*(w = word) = '\0';                        // and reinitialise to empty word
		}
	}
}

// Main function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray* const pistr = prhs[0];                              // pointer to input char array (strings)
	const size_t islen  = mxGetM(pistr);                               // length of input strings
	const size_t nistr  = mxGetN(pistr);                               // number of input strings
mexPrintf("islen = %u, nistr = %u\n",islen,nistr);
	const size_t wmaxlen = ceil(sqrt(2*islen));                        // maximum word length (worst-case scenario!)
	GHashTable* const dict = g_hash_table_new(g_str_hash,g_str_equal); // dictionary implemented as a hash set
	char* const word = malloc(wmaxlen+1);                              // allocate working word cstring buffer
//	for (size_t i=0; i<nistr;++i) {                                    // for each input string in array
//		const char* const istr = mxArrayToString(pistr[i]);       // input cstring buffer
//		LZc(dict,istr,word);                                           // LZ algorithm: build the dictionary
//	}
	free(word);                                                        // deallocate word buffer
	const size_t lzc = g_hash_table_size(dict);                        // LZ complexity = size of dictionary
	plhs[0] = mxCreateDoubleScalar(lzc);                               // output LZ complexity
	if (nlhs > 1) plhs[1] = hash_table_to_cvec(dict);                  // optionally output dictionary as cell vector of strings
	free_hash_table(dict);                                             // clean up dictionary
}
