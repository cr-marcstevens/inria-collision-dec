#!/usr/bin/env python

## \file gen_mmt.py
#Generate code of the sub_isd module using May, Meurer & Thomae algorithm
#
#Parameters p is set as first parameter of the command line. It must multiple of 4.

import os
import sys

def repeat(format, start, end, sep):
	return sep.join([format % (i) for i in range(start, end + 1)])

try:
	ncols=int(sys.argv[1])
except:
	sys.stderr.write('usage : ' +sys.argv[0]+' ncols\n')
	sys.exit(1)

if(ncols%4 != 0):
	sys.stderr.write('p must be multiple of 4')
	sys.exit(2)

output = ""

output += """/* FILE GENERATED BY %s */\n\n""" % os.path.basename(sys.argv[0])

output += """#include <stdio.h>
#include <stdlib.h>
#include "libisd.h"
#include "sub_isd.h"
#include "final_test.h"
#include "measure.h"
#include "waht.h"
#include "iaht.h"
#include "nocolht.h"
"""

output += "#define P %d" % (ncols)

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
	printf("p : %d\\n", params->p);
	printf("eff_word_len : %ld\\n", min(params->r, word_len));
	printf("threshold : %d\\n", params->weight_threshold);
}

void sub_isd_init(isd_params* params, word* local_L, word* local_synds, unsigned int local_N, sw_list** local_h, ranctx* state) {
	params->p = P;
	(void) state;
	print_parameters(params);
	L = local_L;
	N = local_N;
	syndsprime = local_synds;
	h = local_h;

	n = params->n;
	r = params->r;
	k = params->k;
	l = params->l;
	l2 = params->l2;
	w = params->w;

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

	unsigned int i1, i2, i3, i4;
"""
output += "	unsigned int " + repeat("c%d", 1, ncols, ", ") + ';'

output += """
	iaht_list p1, p2, p3, p4;
	unsigned int weight_on_one_word;
	int final_weight;
	word a, aprime, x, old_x, old_sum;
	const word lMSBmask = (1UL << (word_len - l)) - 1; // word with its l MSB to 0
	const int l2shift = word_len - l2;
	const int ll2shift = word_len - (l-l2); 
	word index;
	word value;
	nocolht_elt L12_elt;
	nocolht_elt L34_elt;
	word tmp_value;
	word* E1;
	word* E2;
	word* E3;
	word* E4;
	
"""
output += """/* Intermediate sums */\n"""
for i in range(1, ncols/4+1):
	output += "	word sum" + "".join([str(i) for i in range(1, i+1)]) + ";\n"
#for i in range(1, ncols/4+1):
#	output += "	word sumS" + "".join([str(i) for i in range(ncols/2+1,ncols/2+i+1)]) + ";\n"

output += """
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
output += "		sum1 = L[c1];\n"
for i in range(2,1+ncols/4):
	output += "	for (c%d = %d; c%d < c%d; ++c%d) {\n" % (i, ncols/4-i, i, i-1, i)
	output += "		sum" + "".join([str(j) for j in range(1,i+1)]) + " = sum"  + "".join([str(j) for j in range(1,i)]) + " ^ L[c%d];\n" % (i)
output += "		value = sum" + "".join([str(j) for j in range(1,ncols/4+1)]) +";"

output += """
		index = value >> l2shift;
		waht_store(L1, index, value);
"""
output += "		iaht_store(L1_indices, index, P/4, "+ repeat("c%d", 1, ncols/4, ", ") +");\n	"

for i in range(ncols/4):
	output += "}"
output += "\n\n"

	# L2, L4
output += "	for (c1 = (L_len+1)/2 + %d; c1 < L_len; ++c1) {\n" % (ncols/4-1)
output += "		sum1 = L[c1];\n"
for i in range(2,1+ncols/4):
	output += "	for (c%d = (L_len+1)/2 + %d; c%d < c%d; ++c%d) {\n" % (i, ncols/4-i, i, i-1, i)
	output += "		sum" + "".join([str(j) for j in range(1,i+1)]) + " = sum"  + "".join([str(j) for j in range(1,i)]) + " ^ L[c%d];\n" % (i)
output += "		value = sum" + "".join([str(j) for j in range(1,ncols/4+1)]) +";"
output += """
		index = value >> l2shift;
		waht_store(L2, index, value);
"""
output += "		iaht_store(L2_indices, index, P/4, "+ repeat("c%d", 1, ncols/4, ', ') +");\n	"
for i in range(ncols/4):
	output += "}"
output += "\n\n"

output += """
	for (a = 0; a < (1UL << l2); ++a) {
		aprime = a ^ (synd >> l2shift);
		nocolht_reset(L12, L12_size);

		for (x = 0; x < (1UL << l2); ++x) {

			for(E1 = waht_get(L1, x); E1 != NULL; E1 = waht_next(L1, x, E1)) {
				for(E2 = waht_get(L2, a ^ x); E2 != NULL; E2 = waht_next(L2, a ^ x, E2)) {
					L12_elt = *E1 ^ *E2;
					// the index is the l-l2 bits following the l2 most significant ones
					index = (L12_elt << l2) >> ll2shift;

					// we build the word L12_elt such that it has x on the l2 most significants bits, zero on the l-l2 next bits (wasted space that may be used) and the sum on the rest
					// |   x    | 0 | u+v |
					//  <--l2-->
					//  <----l----->
					//  <----word_len---->

					L12_elt &= lMSBmask; // clear the l MSB
					L12_elt ^= (x << l2shift);  // inject x
					nocolht_store(L12, index, L12_elt);
				}
			}
		}


		for (x = 0; x < (1UL << l2); ++x) {
			for(E3 = waht_get(L3, x); E3 != NULL; E3 = waht_next(L3, x, E3)) {
				tmp_value = *E3 ^ synd;
				for(E4 = waht_get(L4, aprime ^ x); E4 != NULL; E4 = waht_next(L4, aprime ^ x, E4)) {
					L34_elt = tmp_value ^ *E4;
					index = (L34_elt << l2) >> ll2shift;
					L12_elt = nocolht_get(L12, index); // we stored only one element per index in L12
					if (L12_elt != 0) {
						incr_collision_counter();

						value = L12_elt ^ L34_elt;
						value &= lMSBmask; // clear the l MSB; the L12 elt contained the coresponding x in this place (see building of L12)

						weight_on_one_word = isd_weight(value);
						if(weight_on_one_word < threshold) {
							bday_cycle_stopwatch_stop();
							incr_final_test_counter();
							final_test_cycle_stopwatch_start();

							old_x = L12_elt >> l2shift;
							old_sum = L12_elt & lMSBmask;
							/* We found a collision, now we look into the L*_indices to find the corresponding columns numbers */
							for (i1 = 1, p1 = L1_indices[old_x]+1; i1 < 1 + L1_indices[old_x][0]; ++i1, p1 += P/4) {
								for (i2 = 1, p2 = L2_indices[a ^ old_x]+1; i2 < 1 + L2_indices[a ^ old_x][0]; ++i2, p2 += P/4) {
"""
for i in range(1, 1+ncols/2):
	output += """
									c%d = p%d[%d];""" % (i, 1+(i-1)/(ncols/4), (i-1)%(ncols/4))
output += """
									if(((""" + repeat("L[c%d]", 1, ncols/2, " ^ ") + ") & lMSBmask) == old_sum) {\n"
output += """
										for (i3 = 1, p3 = L3_indices[x]+1; i3 < 1 + L3_indices[x][0]; ++i3, p3 += P/4) {
											for (i4 = 1, p4 = L4_indices[aprime ^ x]+1; i4 < 1 + L4_indices[aprime ^ x][0]; ++i4, p4 += P/4) {"""
for i in range(1+ncols/2, 1+ncols):
	output += """
												c%d = p%d[%d];""" % (i, 1+(i-1)/(ncols/4), (i-1)%(ncols/4))
output += """
												if(((""" + repeat("L[c%d]", 1+ncols/2, ncols, " ^ ") + " ^ synd)) == L34_elt) {\n"
output += """
													final_weight = final_test(0, weight_on_one_word, P, """ + repeat("c%d", 1, ncols, ', ') +");"
output += """
													if (final_weight != -1) {"""
output += """
														sw_list_append(h, sw_filled_new(0, final_weight, P, """ + repeat("c%d", 1, ncols, ', ') +"));"
output += """
													}
												}
											}
										}
									}
								}
							}
							final_test_cycle_stopwatch_stop();
							bday_cycle_stopwatch_start();
						}
					}
				}
			}
		}
	}
	sw_list_uniq(h);
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
