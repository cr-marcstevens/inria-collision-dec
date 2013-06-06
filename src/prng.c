#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "prng.h"

uint64_t ranval( ranctx *x ) {
	uint64_t e = x->a - rot(x->b, 7);
	x->a = x->b ^ rot(x->c, 13);
	x->b = x->c + rot(x->d, 37);
	x->c = x->d + e;
	x->d = e + x->a;
	return x->d;
}

void raninit( ranctx *x, uint64_t seed ) {
	uint64_t i;
	x->a = 0xf1ea5eed, x->b = x->c = x->d = seed;
	for (i=0; i<20; ++i) {
		(void)ranval(x);
	}
}

uint64_t random_seed() {
	uint64_t seed;
	FILE* f = fopen("/dev/urandom", "r");
	if (f) {
		if(fread(&seed, sizeof(seed), 1, f) != 1) {
			fprintf(stderr, "Could not read urandom\n");
			seed = time(NULL);
		}
		fclose(f);
	}
	else {
		fprintf(stderr, "Could not open urandom\n");
		seed = time(NULL);
	}
	return seed;
}
