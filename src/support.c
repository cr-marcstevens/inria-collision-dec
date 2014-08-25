#include <stdlib.h>
#include "support.h"

/*
                              L_len
     <----------------------------------------------------->
L = [      a      |      b      |      c      |      d      ]
     <-----------> <-----------> <-----------> <----------->
         L_len/4    (L_len+1)/4   (L_len+2)/4   (L_len+3)/4

half0        = [a|b]
half1        = [c|d]
quarter0and2 = [a|c]
quarter1and3 = [b|d]
 */

#define A (L_len/4)
#define B ((L_len+1)/4)
#define C ((L_len+2)/4)
#define D ((L_len+3)/4)

/* First half */
void prepare_half0(word** Lr, unsigned int* Lr_len, unsigned int L_len) {
	*Lr_len = A + B;
	*Lr = (word*) malloc(*Lr_len * sizeof(word));
}

void build_half0(word* Lr, word* L, unsigned int L_len) {
	memcpy(Lr, L, (A + B) * sizeof(word));
}

ci_t inv_half0(ci_t i, unsigned int L_len) {
	(void) L_len;
	return i;
}

/* Second half */
void prepare_half1(word** Lr, unsigned int* Lr_len, unsigned int L_len) {
	*Lr_len = C + D;
	*Lr = (word*) malloc(*Lr_len * sizeof(word));
}

void build_half1(word* Lr, word* L, unsigned int L_len) {
	memcpy(Lr, L + A + B, (C + D) * sizeof(word));
}

ci_t inv_half1(ci_t i, unsigned int L_len) {
	return i + A + B;
}

/* First and third quarter */
void prepare_quarter0and2(word** Lr, unsigned int* Lr_len, unsigned int L_len) {
	*Lr_len = A + C;
	*Lr = (word*) malloc(*Lr_len * sizeof(word));
}

void build_quarter0and2(word* Lr, word* L, unsigned int L_len) {
	memcpy(Lr, L, A * sizeof(word));
	memcpy(Lr + A, L + A + B, C * sizeof(word));
}

ci_t inv_quarter0and2(ci_t i, unsigned int L_len) {
	if (i < A) {
		return i;
	}
	else {
		return i + B;
	}
}

/* Second and fourth quarter */
void prepare_quarter1and3(word** Lr, unsigned int* Lr_len, unsigned int L_len) {
	*Lr_len = B + D;
	*Lr = (word*) malloc(*Lr_len * sizeof(word));
}

void build_quarter1and3(word* Lr, word* L, unsigned int L_len) {
	memcpy(Lr, L + A, B * sizeof(word));
	memcpy(Lr + B, L + A + B + C, D * sizeof(word));
}

ci_t inv_quarter1and3(ci_t i, unsigned int L_len) {
	if (i < B) {
		return i + A;
	}
	else {
		return i + A + C;
	}
}
