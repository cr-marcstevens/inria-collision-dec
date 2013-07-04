#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libisd.h"
 
double nCr(int n, int r) {
	double x, y, i;
	if(r > n){
		return 0;
	}
	if (r > n/2){
		r = n-r;
	}
	x = 1;
	y = 1;
	i = n-r+1;
	while (i <= n){
		x = (x*i)/y;
		y += 1;
		i += 1;
	}
	return x;
}

void print_bin(word x) {
	unsigned int i;
	for (i = 0; i < word_len; ++i) {
		printf("%lu", (x >> (word_len-i-1)) & 1UL);
	}
}

/* Used for debugging */
void generate_id_permutation(unsigned int* perm, unsigned int* perm_inv, unsigned int n, ranctx* state) {
	(void) state; // remove warning, useful only to keep same prototype as generate_permutation
	unsigned int i;
	fprintf(stderr, "WARNING : NOT PERMUTING\n");
	for (i = 0; i < n; ++i) {
		perm[i] = i;
		perm_inv[i] = i;
	}
}

void generate_permutation(unsigned int* perm, unsigned int* perm_inv, unsigned int n, ranctx* state) {
	unsigned int i, j, tmp;
	for (i = 0; i < n; ++i) {
		perm[i] = i;
	}
	for (i = n-1; i >= 1; --i) {
		j = random_range(state, i+1);
		tmp = perm[j]; perm[j] = perm[i]; perm[i] = tmp;
	}
	if (perm_inv != NULL) {
		for (i = 0; i < n; ++i) {
			perm_inv[perm[i]] = i;
		}
	}
}

void apply_permutation(mzd_t* dst, mzd_t* src, unsigned int* perm, unsigned int n) {
	unsigned int i;
	for (i = 0; i < n; ++i) {
		mzd_copy_row (dst, perm[i], src, i);
	}
}

unsigned int isd_parity(word w) {
	return __builtin_parityll(w);
}

unsigned int isd_weight(word w) {
	return __builtin_popcountll(w);
}
