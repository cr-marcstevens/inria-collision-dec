#include <stdio.h>
#include <stdlib.h>
#include "libisd.h"
#include "process_solutions.h"
#include <m4ri/m4ri.h>

/** function called called on each solution an iteration could find
 * It takes a pointer to sparse word of weight p as an input, 
 * then
 *   multiplies it by BT, 
 *   adds the syndrome, 
 *   obtains the coresponding vector of weight w-p, 
 *   write on stdout the w columns solutions of the problem. 
 */

void process_solution(sw_list* eprime, unsigned int w, unsigned int l, mzd_t* BT, word** syndzero, mzd_t* U, unsigned int* perm_inv, unsigned long long nb_iter) {
	unsigned int i, j;
	unsigned int sols_idx;
	unsigned int r = BT->ncols;
	word chunk;
	sw* complete_sol = sw_new(w);
	word* eprimeBT = (word*) malloc(bit_in_words(r) * sizeof(word));

	// compute syndzero + e'*BT
	for (i = 0; i < bit_in_words(r); ++i) {
		eprimeBT[i] = syndzero[eprime->synd_idx][i];
	}

	for (i = 0; i < eprime->weight; ++i) {
		for (j = 0; j < bit_in_words(r); ++j) {
			eprimeBT[j] ^= BT->rows[eprime->pos[i]][j];
		}
	}

	// add the p positions to the complete_sol vector
	sols_idx = 0;
	for (i = 0; i < eprime->weight; ++i) {
		complete_sol->pos[sols_idx] = perm_inv[r-l+eprime->pos[i]];
		++sols_idx;
	}

	// build each bits of (syndzero + e'*BT)*U one by one, each non null bit corespond to a position of the complete_sol.
	for (i = 0; i < r; ++i) {
		chunk = 0;
		for (j = 0; j < bit_in_words(r); ++j) {
			chunk ^= eprimeBT[j] & U->rows[i][j];
		}
		if (isd_parity(chunk)) {
			complete_sol->pos[sols_idx] = perm_inv[i];
			++sols_idx;
		}
	}
	complete_sol->synd_idx = eprime->synd_idx;
	complete_sol->weight = sols_idx; // the algorithm may have found a word with weight lower than w-p

	sw_sort(complete_sol);

	printf("%12lld ", nb_iter);
	sw_print(complete_sol);

	free(eprimeBT);
	sw_free(complete_sol);
}
