#include "mex.h"
#include <string.h>

/* n = LZW_mex(s)

   WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING AT ALL! CALL IT RIGHT!

   Return indices of first occurrences of elements of row vector a in row
   vector b, in row vector idx. Assumes elements are integers, since double
   float comparisons are exact. Returned indices are 1-offset, or zero if
   element of a not found. (Dumb linear search satisfices.)

	mex -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Wall -Werror -O3" C/LZW_mex.c -outdir mex

   To compile, see mvgc_makemex.m in the utils subfolder */

// Minimal singly-linked list implementation

typedef struct sll
{
	char* pdata;
	struct sll* pnext;
} sll_t;

sll_t* const sll_add_data(sll_t* psll, const char* const pdata, const size_t n) // IMPORTANT: size n MUST INCLUDE terminating null byte
{
	sll_t* posll = psll;
	psll = malloc(sizeof(sll_t));
	psll->pnext = posll;
	psll->pdata = malloc(n);
	memcpy(psll->pdata,pdata,n); // since we have length, memcpy is more efficient than strcpy
	return psll;
}

void sll_free(sll_t* psll)
{
	while (psll != NULL) {
		sll_t* pnsll = psll->pnext;
		free(psll->pdata);
		free(psll);
		psll = pnsll;
	}
}

void sll_list(const sll_t* psll)
{
	while (psll != NULL) {
		mexPrintf("'%s'\n",psll->pdata);
		psll = psll->pnext;
	}
}

bool sll_find(const sll_t* psll,const char* const str)
{
	while (psll != NULL) {
		if (strcmp(str,psll->pdata) == 0) {
			return true;
		}
		psll = psll->pnext;
	}
	return false;
}

size_t sll_size(const sll_t* psll)
{
	size_t n = 0;
	while (psll != NULL) {
		++n;
		psll = psll->pnext;
	}
	return n;
}

// Main function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray* const pistr = prhs[0]; // pointer to input string

	const char* const s = mxArrayToString(pistr); // input string

	char* const word = malloc(mxGetM(pistr)*mxGetN(pistr)+1); // length of input string: overkill, but safe!

	sll_t* dict = NULL; // the dictionary

	word[0] = *s;
	size_t n = 1; // length of current word
	for (const char* c = s+1; *c != '\0'; ++c) {
		word[n] = *c;
		word[++n] = '\0';
		if (!sll_find(dict,word)) { // if word not found, add to dictionary and reset to current character
			dict = sll_add_data(dict,word,n+1); // need to supply length+1 to sll_add_data
			word[0] = *c;
			word[n = 1] = '\0';
		}
	}

	plhs[0] = mxCreateDoubleScalar(sll_size(dict)); // output = size of dictionary

	sll_list(dict);

	sll_free(dict);

	free(word);
}
