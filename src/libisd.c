#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libisd.h"
 
double nCr(int n, int r) {
	double x, y, i;
	if (r < 0) {
		return 0;
	}
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

int word_cmp (const void * a, const void * b) {
	if (*(word*)a < *(word*)b) return -1;
	if (*(word*)a > *(word*)b) return 1;
	return 0;
}

int pscal(word* a, word* b, int nb_bits) {
	unsigned int j;
	word chunk = 0;
	for (j = 0; j < bit_in_words(nb_bits); ++j) {
		chunk ^= a[j] & b[j];
	}
	return isd_parity(chunk);
}

unsigned int isd_parity(word w) {
	return __builtin_parityll(w);
}

unsigned int isd_weight(word w) {
	return __builtin_popcountll(w);
}
