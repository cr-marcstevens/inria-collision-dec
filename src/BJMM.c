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
#include "BJMMtools.h"

static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, l2, l3, p, e1, e2, w, L_len, threshold, csize;
static unsigned int p2;
static unsigned int p1;
static int shift;
static sw_list** h;

void sub_isd_init(word* simple_HprimemodT, unsigned int local_N, word* local_syndsprime, unsigned int local_n, unsigned int local_r,unsigned int local_l, unsigned int local_l2, unsigned int local_l3, unsigned int local_p, unsigned int local_e1, unsigned int local_e2, unsigned int local_w, unsigned int local_threshold,unsigned int local_csize, sw_list** local_h) {
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
	csize= local_csize;

	k = n-r;

	L_len = k+l;

	threshold = local_threshold;
	shift = min(r, word_len) - l;

	p1= p/2 + e1;
	p2 = p1/2 + e2;
}


void sub_isd() {

	/*
	 * construction des listes  a partir de E1, E2, E3, E4
	 */

	short* indice = malloc((p2/2)*sizeof(short));
	word** sums = malloc((p2/2)*csize*sizeof(word));


	//TODO

	free(indice);
	free(sums);
}

void sub_isd_report(unsigned long long cycles_per_iter) {
	(void) cycles_per_iter;
}

void sub_isd_free() {

}
