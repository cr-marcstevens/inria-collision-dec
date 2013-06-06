#ifndef PRNG_H
#define PRNG_H

#include <stdint.h>

/* random number generator with a 64 bits seed : from http://burtleburtle.net/bob/rand/smallprng.html */
typedef struct ranctx {
	uint64_t a;
	uint64_t b;
	uint64_t c;
	uint64_t d;
} ranctx;

#define rot(x,k) (((x)<<(k))|((x)>>(64-(k))))
uint64_t ranval( ranctx *x );
void raninit( ranctx *x, uint64_t seed );
uint64_t random_seed();
//#define random_range(RANGE) (random() / (RAND_MAX / RANGE + 1))
#define random_range(state, RANGE) (ranval(state) % (RANGE))

#endif
