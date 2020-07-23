#include "mex.h"
#include <string.h>
#include <gmodule.h>

/* n = LZW_mex(s)

   WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING AT ALL! CALL IT RIGHT!

   Return indices of first occurrences of elements of row vector a in row
   vector b, in row vector idx. Assumes elements are integers, since double
   float comparisons are exact. Returned indices are 1-offset, or zero if
   element of a not found. (Dumb linear search satisfices.)

	mex -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Wall -Werror -O3" C/LZW_mex.c -outdir mex

   To compile, see mvgc_makemex.m in the utils subfolder */

#define NUL ('\0')

//GHashTable* htab = g_hash_table_new(g_str_hash,GEqualFunc key_equal_func);


// Minimal singly-linked list implementation of dictionary (linear search!)

void printstr(const char* const str,const size_t n) // non-null-terminated!
{
	char* const cstr = malloc(n+1);
	memcpy(cstr,str,n);
	cstr[n] = NUL;
	mexPrintf("'%s'",cstr);
	free(cstr);
}

typedef struct sll
{
	char* str;
	size_t slen;
	struct sll* next;
} dict_t;

dict_t* const dict_add_data(dict_t* dict, const char* const str, const size_t n)
{
	dict_t* odict = dict;
	dict = malloc(sizeof(dict_t));
	dict->str = malloc(n);
	memcpy(dict->str,str,n);
	dict->slen = n;
	dict->next = odict;
	return dict;
}

void dict_free(dict_t* dict)
{
	while (dict != NULL) {
		dict_t* ndict = dict->next;
		free(dict->str);
		free(dict);
		dict = ndict;
	}
}

void dict_list(const dict_t* dict)
{
	while (dict != NULL) {
		printstr(dict->str,dict->slen);
		mexPrintf("\n");
		dict = dict->next;
	}
}

bool dict_find(const dict_t* dict, const char* const str, const size_t n)
{
	while (dict != NULL) {
		if (dict->slen == n) {
			if (memcmp(str,dict->str,n) == 0) {
				return true;
			}
		}
		dict = dict->next;
	}
	return false;
}

size_t dict_size(const dict_t* dict)
{
	size_t N = 0;
	while (dict != NULL) {
		++N;
		dict = dict->next;
	}
	return N;
}

// Main function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const char* s = mxArrayToString(prhs[0]); // input string buffer

	dict_t* dict = NULL; // the dictionary

	const char* w = s;
	for (++s;*s;++s) { // terminating condition equiv to *s == NUL !!!
		const size_t n = s-w+1; // length of current word
		if (!dict_find(dict,w,n)) { // if word not found, add to dictionary and reset to current character
			dict = dict_add_data(dict,w,n);
			w = s;
		}
	}

	const size_t N = dict_size(dict);

	plhs[0] = mxCreateDoubleScalar(N); // output = size of dictionary

	// dict_list(dict);

	if (nlhs > 1) { // output dictionary as cell array
		mxArray* cmat = mxCreateCellMatrix(N,1);
		const dict_t* pdict = dict; // the dictionary
		size_t i = N;
		while (pdict != NULL) {
			const size_t n = pdict->slen;
			char* const cstr = malloc(n+1);
			memcpy(cstr,pdict->str,n);
			cstr[n] = NUL;
			mxSetCell(cmat,--i,mxCreateString(cstr));
			free(cstr);
			pdict = pdict->next;
		}
		plhs[1] = cmat;
	}

	dict_free(dict);
}

#undef NUL
