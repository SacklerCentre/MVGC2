#include "mex.h"
#include <string.h>
#include <gmodule.h>

/*
	Calculate LZW complexity, and (optionally) return the dictionary. Call in Matlab as (see LWZ.m)

		[n,dict] = LZW(s,true);

	NOTE: This is C99, and needs Glib "gmodule" for GHashTable.

	To compile on Linux with GCC (static linkage for portability):

		mex -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Wall -Werror -O3 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -static" -lglib-2.0 LZW_mex.c

	This should also work on Windows with MinGW-w64 (https://mingw-w64.org/doku.php),
	(note: MS Visual Studio doesn't do C99), and probably fine on Mac, but not tested.

		*** WARNING: out-of-memory not checked (TODO)! ***

	Benchmark: for a 1,000,000 character binary input string:

		Matlab  (using containers.Map) : 7.31 seconds.
		LZW_mex (using GHashTable)     : 0.11 seconds

*/

/*
void print_hash_table(GHashTable* dict)
{
	GHashTableIter iter;
	gpointer key;
	g_hash_table_iter_init(&iter,dict);
	while (g_hash_table_iter_next(&iter,&key,NULL)) mexPrintf("\t'%s'\n",key);
}
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

const size_t LZW(const char* const istr, GHashTable* const dict, const size_t wmaxlen)
{
	char* const word = malloc(wmaxlen+1);              // allocate working word string buffer
	const char* c = istr;                              // initialise current input character to start of input string
	size_t n;                                          // word length
	word[0] = *c; word[n=1] = '\0';                    // initialise word to first input character
	for (++c;*c;++c) {                                 // terminating condition equiv to *c == '\0' !!!
		word[n] = *c;                                  // append next input character to word
		if (++n > wmaxlen) {                           // buffer overflow alert!
			free(word);                                // deallocate word buffer
			return -1;                                 // error condition (note unsigned wraps!)
		}
		word[n] = '\0';
		if (g_hash_table_lookup(dict,word) == NULL) {  // if word not found,
			g_hash_table_add(dict,strdup(word));       // add to dictionary
			word[0] = *c; word[n=1] = '\0';            // and reset word to current input character;
		}
	}
	if (n > 1) {                                       // is there a word "left over"?
		if (g_hash_table_lookup(dict,word) == NULL) {  // if word not in dictionary,
			g_hash_table_add(dict,strdup(word));       // add to dictionary
		}
	}
	free(word);                                        // deallocate word buffer
	return g_hash_table_size(dict);                    // return LZW complexity (size of dictionary)
}


// Main function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const char* const istr = mxArrayToString(prhs[0]);                 // input string buffer

	GHashTable* const dict = g_hash_table_new(g_str_hash,g_str_equal); // dictionary implemented as a hash set

	const size_t wmaxlen = mxGetScalar(prhs[1]);                       // max. word length

	const size_t lzwc = LZW(istr,dict,wmaxlen);                        // LZW algorithm: build the dictionary

	if (lzwc == -1) { // error (word length too small)
		plhs[0] = mxCreateDoubleScalar(-1);                            // output -1
		if (nlhs > 1) plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);      // optionally output empty matrix
	}
	else {
		plhs[0] = mxCreateDoubleScalar(lzwc);                          // output LZW complexity
		if (nlhs > 1) plhs[1] = hash_table_to_cvec(dict);              // optionally output dictionary as cell vector of strings
	}

	free_hash_table(dict); // clean up dictionary
}
