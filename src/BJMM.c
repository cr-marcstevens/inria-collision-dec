/*
 * BJMM.c
 *
 *  Created on: 25 juin 2013
 *      Author: Mathieu ARIA
 */

//TODO

#include <stdlib.h>
#include <time.h>
#include "sub_isd.h"
#include "m4ri/m4ri.h"
#include "libisd.h"
#include "sparse_words_list.h"
#include "final_test.h"
#include "measure.h"
#include "cpucycles/cpucycles.h"

static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, l2, l3, p, e1, e2, w, L_len, threshold;
static int shift;
static sw_list** h;

void sub_isd_init(word* simple_HprimemodT, unsigned int local_N, word* local_syndsprime, unsigned int local_n, unsigned int local_r,unsigned int local_l, unsigned int local_l2, unsigned int local_l3, unsigned int local_p, unsigned int local_e1, unsigned int local_e2, unsigned int local_w, unsigned int local_threshold, sw_list** local_h) {
	L = simple_HprimemodT;
	N = local_N;
	syndsprime = local_syndsprime;
	n = local_n;
	r = local_r;
	l = local_l;
	l2 = local_l2; //  named r2 in the paper
	l3 = local_l3; //  named r1 in the paper
	e1 = local_e1; //  filter 1 parameter
	e2 = local_e2; //  filter 2 parameter
	p = local_p;
	w = local_w;
	h = local_h;

	k = n-r;

	L_len = k+l;

	threshold = local_threshold;
	shift = min(r, word_len) - l;
}


void sub_isd() {
	word synd = syndsprime[0]; //DOOM not implemented
	(void) synd;
	//TODO
}

void sub_isd_report(unsigned long long cycles_per_iter) {
	(void) cycles_per_iter;
}

void sub_isd_free() {
	//TODO
}
