#include <stdio.h>
#include <stdlib.h>
#include "libisd.h"
#include "process_solutions.h"
#include "m4ri/m4ri.h"


/* function called for each found solution. It takes a sparse word of weight p as an input, multiply it by BT, adds the syndrome, obtains the coresponding vector of weight w-p and write on stdout the w columns solutions of the problem. */

void process_solutions_on_the_fly(sw_list** h, unsigned int w, unsigned int l, mzd_t* BT, word** syndzero, mzd_t* U, unsigned int* perm_inv, unsigned long long nb_iter) {
	unsigned int i, j;
	unsigned int sols_idx;
	unsigned int r = BT->ncols;
	word chunk;
	sw_list* solution = sw_list_new(NULL, w);
	sw_list* eprime = *h;
	word* eprimeBT = (word*) malloc(bit_in_words(r) * sizeof(word));

	while(eprime) {
		// compute syndzero + e'*BT
		for (i = 0; i < bit_in_words(r); ++i) {
			eprimeBT[i] = syndzero[eprime->synd_idx][i];
		}

		for (i = 0; i < eprime->weight; ++i) {
			for (j = 0; j < bit_in_words(r); ++j) {
				eprimeBT[j] ^= BT->rows[eprime->pos[i]][j];
			}
		}

		// add the p positions to the solution vector
		sols_idx = 0;
		for (i = 0; i < eprime->weight; ++i) {
			solution->pos[sols_idx] = perm_inv[r-l+eprime->pos[i]];
			++sols_idx;
		}

		// build each bits of (syndzero + e'*BT)*U one by one, each non null bit corespond to a position of the solution.
		for (i = 0; i < r; ++i) {
			chunk = 0;
			for (j = 0; j < bit_in_words(r); ++j) {
				chunk ^= eprimeBT[j] & U->rows[i][j];
			}
			if (isd_parity(chunk)) {
				solution->pos[sols_idx] = perm_inv[i];
				++sols_idx;
			}
		}
		solution->synd_idx = eprime->synd_idx;
		solution->weight = sols_idx; // the algorithm may have found a word with weight lower than w-p

		sw_list_sort(solution);
		printf("%12lld ", nb_iter);
		
		sw_list_print(solution);

		eprime = eprime->next;
	}
	sw_list_free(solution);
	sw_list_free(*h);
	*h = NULL;
}


void process_solutions_at_end(sw_list** h, unsigned int w, unsigned int l, mzd_t* BT, word** syndzero, mzd_t* U, unsigned int* perm_inv) {
	(void) h;
	(void) w;
	(void) l;
	(void) BT;
	(void) syndzero;
	(void) U;
	(void) perm_inv;
	return;
}
