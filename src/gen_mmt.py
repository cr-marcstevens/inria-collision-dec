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
void sub_isd_init(word* simple_HprimemodT, unsigned int local_N, word* local_syndsprime, unsigned int local_n, unsigned int local_r,unsigned int local_l, unsigned int local_l2, unsigned int local_l3, unsigned int local_p, unsigned int local_e1, unsigned int local_e2, unsigned int local_w, unsigned int local_threshold, unsigned int local_csize, sw_list** local_h) {
	L = simple_HprimemodT;
	N = local_N;
	syndsprime = local_syndsprime;
	n = local_n;
	r = local_r;
	l = local_l;
	l2 = local_l2;
	(void) local_l3;
	(void) local_e1;
	(void) local_e2;
	(void) local_p;
	(void) local_csize;
	w = local_w;
	h = local_h;

	k = n-r;

	L_len = k+l;

	threshold = local_threshold;

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
	word lmask = (1UL << l) - 1;
	word l2mask = (1UL << l2) - 1;
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
		index = value & l2mask;
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
		index = value & l2mask;
		waht_store(L2, index, value);
"""
output += "		iaht_store(L2_indices, index, p/4, "+ repeat("c%d", ncols/4, ', ') +");\n	"
for i in range(ncols/4):
	output += "}"
output += "\n\n"

output += """
	for (a = 0; a < (1UL << l2); ++a) {
		aprime = (a ^ synd) & l2mask;
		nocolht_reset(L12, L12_size);

		for (x = 0; x < (1UL << l2); ++x) {

			waht_list E1 = waht_get(L1, x);
			if (E1 != NULL) {
				waht_list E2 = waht_get(L2, a ^ x);
				if (E2 != NULL) {
					for (i = 1; i < E1[0]; ++i) {
						for (j = 1; j < E2[0]; ++j) {
							value = E1[i] ^ E2[j];
							// the index is the l first bits of the sum minus the l2 first
							index = (value & lmask) >> l2;

							// we build the word value such that it has x on the l2 first bits, zero on the l-l2 next bits (wasted space that may be used) and the sum on the rest
							// |   x    | 0 | u+v |
							//  <--l2-->
							//  <----l----->

							value &= ~lmask;
							value ^= x;
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
							index = (value & lmask) >> l2;
							nocolht_elt L12_elt = nocolht_get(L12, index); // we stored only one element per index in L12
							if (L12_elt != 0) {
								incr_nb_collision_counter();

								value ^= L12_elt;
								value >>= l; // the l2 first bits of L12_elt are unknown (but the l first bits of the sum is 0 at this point of the algorithm)

								if(isd_weight(value) <= threshold) {
									old_x = L12_elt & l2mask;

									/* Ok we found a collision, now we look into the L*_indices to find the corresponding columns numbers */
									for (i1 = 1, p1 = L1_indices[old_x]+1; i1 < 1 + L1_indices[old_x][0]; ++i1, p1 += p/4) {
										for (i2 = 1, p2 = L2_indices[a ^ old_x]+1; i2 < 1 + L2_indices[a ^ old_x][0]; ++i2, p2 += p/4) {
											for (i3 = 1, p3 = L3_indices[x]+1; i3 < 1 + L3_indices[x][0]; ++i3, p3 += p/4) {
												for (i4 = 1, p4 = L4_indices[aprime ^ x]+1; i4 < 1 + L4_indices[aprime ^ x][0]; ++i4, p4 += p/4) {
"""
for i in range(1, 1+ncols):
	output += "													c%d = p%d[%d];\n" % (i, 1+(i-1)/(ncols/4), (i-1)%(ncols/4))
output += "													weight_on_one_word = isd_weight(synd ^ " + repeat("L[c%d]", ncols, "  ^ ") + ");"
output+= """
													if(weight_on_one_word <= threshold) {
														incr_final_test_counter();
														final_test_probe_start();
"""
output += "														final_weight = final_test(0, weight_on_one_word, p, " + repeat("c%d", ncols, ', ') +");"
output += """
														final_test_probe_stop();

														if (final_weight != -1) {
															"""
output += "															*h = sw_list_add(*h, 0, final_weight, p, " + repeat("c%d", ncols, ', ') +");"
output += """
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


void sub_isd_report(unsigned long long cycles_per_iter) {
	(void) cycles_per_iter;
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