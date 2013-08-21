#!/usr/bin/env python

import sys

def repeat(format, nb, sep):
	return sep.join([format % (i) for i in range(1, nb + 1)])

try:
	ncols=int(sys.argv[1])
except:
	sys.stderr.write('usage : ' +sys.argv[0]+' ncols\n')
	sys.exit(1)

if(ncols%4 != 0):
	sys.stderr.write('p must be multiple of 4')
	sys.exit(2)

output = ""

output += """#include <stdlib.h>
#include "libisd.h"
#include "sub_isd.h"
#include "final_test.h"
#include "measure.h"
#include "waht.h"
#include "iaht.h"
#include "nocolht.h"
"""

output += "#define p %d" % (ncols)

output += """
static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, l2, w, L_len, threshold;
static sw_list** h;


static unsigned int L1_size, L2_size, L3_size, L4_size, L12_size, L1_indices_size, L2_indices_size, L3_indices_size, L4_indices_size;
static waht L1 = NULL;
static waht L2 = NULL;
static waht L3 = NULL;
static waht L4 = NULL;
static iaht L1_indices = NULL;
static iaht L2_indices = NULL;
static iaht L3_indices = NULL;
static iaht L4_indices = NULL;
static nocolht L12 = NULL;
"""

output += """
void print_parameters(isd_params* params) {
	printf("n : %d\\n", params->n);
	printf("r : %d\\n", params->r);
	printf("w : %d\\n", params->w);
	printf("l : %d\\n", params->l);
	printf("l2 : %d\\n", params->l2);
	printf("p : %d\\n", p);
	printf("eff_word_len : %ld\\n", min(params->r, word_len));
	printf("threshold : %d\\n", params->weight_threshold);
}

void sub_isd_init(isd_params* params, word* local_L, word* local_synds, unsigned int local_N, sw_list** local_h, ranctx* state) {
	(void) state;
	print_parameters(params);
	L = local_L;
	N = local_N;
	syndsprime = local_synds;
	h = local_h;

	n = params->n;
	r = params->r;
	l = params->l;
	l2 = params->l2;
	w = params->w;

	k = n-r;

	L_len = k+l;

	threshold = params->weight_threshold;

	L1_size = 1UL << l2;
	L2_size = 1UL << l2;
	L3_size = 1UL << l2;
	L4_size = 1UL << l2;
	L1_indices_size = 1UL << l2;
	L2_indices_size = 1UL << l2;
	L3_indices_size = 1UL << l2;
	L4_indices_size = 1UL << l2;
	L12_size = 1ULL<<(l-l2);

	L1 = waht_init(L1_size);
	L2 = waht_init(L2_size);
	L3 = L1;
	L4 = L2;
	//L3 = waht_init(L3_size);
	//L4 = waht_init(L4_size);

	L1_indices = iaht_init(L1_indices_size);
	L2_indices = iaht_init(L2_indices_size);
	L3_indices = L1_indices;
	L4_indices = L2_indices;
	//L3_indices = iaht_init(L3_indices_size);
	//L4_indices = iaht_init(L4_indices_size);

	L12 = nocolht_init(L12_size);
}
"""

output += """
void sub_isd() {
	word synd = syndsprime[0]; //DOOM not implemented

	unsigned int i, j, i1, i2, i3, i4;
"""
output += "	unsigned int " + repeat("c%d", ncols, ", ") + ';'

output += """
	iaht_list p1, p2, p3, p4;
	unsigned int weight_on_one_word;
	int final_weight;
	word a, aprime, x, old_x;
	word lmask = (1UL << (word_len - l)) - 1;
	int l2shift = word_len - l2;
	int ll2shift = word_len - (l-l2); 
	word index;
	word value;

	waht_reset(L1, L1_size);
	waht_reset(L2, L2_size);
	//waht_reset(L3, L3_size);
	//waht_reset(L4, L4_size);
	iaht_reset(L1_indices, L1_indices_size);
	iaht_reset(L2_indices, L2_indices_size);
	//iaht_reset(L3_indices, L3_indices_size);
	//iaht_reset(L4_indices, L4_indices_size);

"""
	# L1, L3
output += "	for (c1 = %d; c1 < L_len/2; ++c1) {\n" % (ncols/4-1)
for i in range(2,1+ncols/4):
	output += "	for (c%d = %d; c%d < c%d; ++c%d) {\n" % (i, ncols/4-i, i, i-1, i)
output += "		value = " + repeat("L[c%d]", ncols/4, " ^ ") + ';'

output += """
		index = value >> l2shift;
		waht_store(L1, index, value);
"""
output += "		iaht_store(L1_indices, index, p/4, "+ repeat("c%d", ncols/4, ", ") +");\n	"

for i in range(ncols/4):
	output += "}"
output += "\n\n"

	# L2, L4
output += "	for (c1 = (L_len+1)/2 + %d; c1 < L_len; ++c1) {\n" % (ncols/4-1)
for i in range(2,1+ncols/4):
	output += "	for (c%d = (L_len+1)/2 + %d; c%d < c%d; ++c%d) {\n" % (i, ncols/4-i, i, i-1, i)
output += "		value = " + repeat("L[c%d]", ncols/4, " ^ ") + ';'
output += """
		index = value >> l2shift;
		waht_store(L2, index, value);
"""
output += "		iaht_store(L2_indices, index, p/4, "+ repeat("c%d", ncols/4, ', ') +");\n	"
for i in range(ncols/4):
	output += "}"
output += "\n\n"

output += """
	for (a = 0; a < (1UL << l2); ++a) {
		aprime = a ^ (synd >> l2shift);
		nocolht_reset(L12, L12_size);

		for (x = 0; x < (1UL << l2); ++x) {

			waht_list E1 = waht_get(L1, x);
			if (E1 != NULL) {
				waht_list E2 = waht_get(L2, a ^ x);
				if (E2 != NULL) {
					for (i = 1; i < E1[0]; ++i) {
						for (j = 1; j < E2[0]; ++j) {
							value = E1[i] ^ E2[j];
							// the index is the l-l2 bits following the l2 most significant ones
							index = (value << l2) >> ll2shift;

							// we build the word value such that it has x on the l2 most significants bits, zero on the l-l2 next bits (wasted space that may be used) and the sum on the rest
							// |   x    | 0 | u+v |
							//  <--l2-->
							//  <----l----->
							//  <----word_len---->

							value &= lmask; // clear the l MSB
							value ^= (x << l2shift);  // inject x
							nocolht_store(L12, index, value);
						}
					}
				}
			}
		}


		for (x = 0; x < (1UL << l2); ++x) {
			waht_list E3 = waht_get(L3, x);
			if (E3 != NULL) {
				waht_list E4 = waht_get(L4, aprime ^ x);
				if (E4 != NULL) {
					for (i = 1; i < 1+E3[0]; ++i) {
						for (j = 1; j < 1+E4[0]; ++j) {
							value = E3[i] ^ E4[j] ^ synd;
							index = (value << l2) >> ll2shift;
							nocolht_elt L12_elt = nocolht_get(L12, index); // we stored only one element per index in L12
							if (L12_elt != 0) {
								incr_nb_collision_counter();

								value ^= L12_elt;
								value = (value << l) >> l; // clear the l MSB; the L12 contained the coresponding x in this place (see building of L12)

								weight_on_one_word = isd_weight(value);
								if(weight_on_one_word <= threshold) {
									incr_final_test_counter();
									final_test_probe_start();
									old_x = L12_elt >> l2shift;

									/* We found a collision, now we look into the L*_indices to find the corresponding columns numbers */
									for (i1 = 1, p1 = L1_indices[old_x]+1; i1 < 1 + L1_indices[old_x][0]; ++i1, p1 += p/4) {
										for (i2 = 1, p2 = L2_indices[a ^ old_x]+1; i2 < 1 + L2_indices[a ^ old_x][0]; ++i2, p2 += p/4) {
											for (i3 = 1, p3 = L3_indices[x]+1; i3 < 1 + L3_indices[x][0]; ++i3, p3 += p/4) {
												for (i4 = 1, p4 = L4_indices[aprime ^ x]+1; i4 < 1 + L4_indices[aprime ^ x][0]; ++i4, p4 += p/4) {"""
for i in range(1, 1+ncols):
	output += """
													c%d = p%d[%d];""" % (i, 1+(i-1)/(ncols/4), (i-1)%(ncols/4))
output += """
													final_weight = final_test(0, weight_on_one_word, p, """ + repeat("c%d", ncols, ', ') +");"
output += """

													if (final_weight != -1) {"""
output += """
														*h = sw_list_add(*h, 0, final_weight, p, """ + repeat("c%d", ncols, ', ') +");"
output += """
													}
												}
											}
										}
									}
									final_test_probe_stop();
								}
							}
						}
					}
				}
			}
		}
	}
}


void sub_isd_report(unsigned long long cycles_per_iter, long long pivot_cost, long long bday_cost, long long final_test_cost) {
	(void) cycles_per_iter;
	(void) pivot_cost;
	(void) bday_cost;
	(void) final_test_cost;
}

void sub_isd_free() {
	waht_free(L1, L1_size);
	waht_free(L2, L2_size);
	//waht_free(L3, L3_size);
	//waht_free(L4, L4_size);
	iaht_free(L1_indices, L1_indices_size);
	iaht_free(L2_indices, L2_indices_size);
	//iaht_free(L3_indices, L3_indices_size);
	//iaht_free(L4_indices, L4_indices_size);
	nocolht_free(L12);
}"""

print output
