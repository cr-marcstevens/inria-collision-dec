//#include <stdio.h>
//#include <stdlib.h>
//#include <stdint.h>
//#include <signal.h>
#include <time.h>
#include "isd.h"
#include "m4ri/m4ri.h"
#include "custom_brilliantrussian.h"
#include "sub_isd.h"
#include "libisd.h"
#include "final_test.h"
#include "process_solutions.h"
#include "sparse_words_list.h"
//#include "io.h"
#include "prng.h"
#include "cpucycles/cpucycles.h"

#include "measure.h"

static int print_status = 0;
static int stop = 0;

void status_handler() {
	print_status = 1;
}

void stop_handler() {
	stop = 1;
}

void scramble(mzd_t* H, unsigned int N, word** synds) {
	int r = H->ncols;
	mzd_t* M = mzd_init(r, r);
	mzd_t* Mcopy = mzd_init(r, r);
	mzd_t* sH = mzd_copy(NULL, H);
	do {
		mzd_randomize(M);
		// mzd_echelonize modifies M so we work on a copy to keep M random
		mzd_copy(Mcopy, M); 
	}	while (mzd_echelonize(Mcopy, 0) != r);
	mzd_mul_m4rm(sH, H, M, 0);
	mzd_copy(H, sH);
	mzd_free(Mcopy);
	mzd_free(sH);

	unsigned int i;
	int j;
	mzd_t* s = mzd_init(1, r);
	mzd_t* ss = mzd_init(1, r);
	for (i = 0; i < N; ++i) {
		for (j = 0; j < s->width; ++j) {
			s->rows[0][j] = synds[i][j];
		}
		mzd_mul_m4rm(ss, s, M, 0);
		for (j = 0; j < s->width; ++j) {
			synds[i][j] = ss->rows[0][j];
		}
	}
	mzd_free(s);
	mzd_free(ss);

	mzd_free(M);
}

// concatenate A and an identity matrix (mzd_concat is slow) into A_I. A_I must be allocated. Does not check for dimensions.
void mzd_concat_id(mzd_t* A_I, mzd_t* A) {
	int i;
	mzd_set_ui(A_I, 0);
	mzd_copy(A_I, A);
	for (i = 0; i < A_I->nrows; ++i) {
		mzd_write_bit(A_I, i, i+A->ncols, 1);
	}
}


int mzd_partial_echelonize(mzd_t* A, int l) {
	if(_mzd_partial_echelonize_m4ri(A, 1, 0, 0,0, A->nrows-l) != A->nrows-l) {
		return 1;
	}
	else {
		return 0;
	}
}

/* This computes an echelonization on the first l columns of A */
/*
int mzd_partial_echelonize(mzd_t* A, int l) {
	int n = A->ncols;
	int r = A->nrows;
	mzd_t *A00;
	mzd_t *A10;
	mzd_t *A01;
	mzd_t *A11;
	mzd_t *A0;
	mzd_t *A1;
	mzp_t *P;
	mzp_t *Q;
	int rank;

	A0 = mzd_init_window(A, 0, 0, r, r-l);
	P = mzp_init(A0->nrows);
	Q = mzp_init(A0->ncols);

	rank = mzd_ple(A0, P, Q, 0);

	if (rank != r-l) { //Echelonisation failed
		mzd_free_window(A0);
		mzp_free(P);
		mzp_free(Q);
		return 1;
	}

	A1 = mzd_init_window(A, 0, r-l, r, n);

	A00 = mzd_init_window(A, 0, 0, r-l, r-l);
	A10 = mzd_init_window(A, r-l, 0, r, r-l);
	A01 = mzd_init_window(A, 0, r-l, r-l, n);
	A11 = mzd_init_window(A, r-l, r-l, r, n);

	mzd_apply_p_left(A1, P);
	_mzd_trsm_lower_left(A00, A01, 0);


	mzd_addmul(A11, A10, A01, 0);

	mzd_apply_p_right_trans_tri(A00, Q);

	// update non-pivot columns

	mzd_trsm_upper_left(A00, A01, 0);

	mzd_set_ui(A00, 0);
	for(rci_t i = 0; i < r-l; ++i)
		mzd_write_bit(A, i, i, 1);

	mzd_free_window(A00);

	// expand E again
	A00 = mzd_init_window(A0, 0, 0, r-l, A0->ncols);
	mzd_apply_p_right(A00, Q);

	// clear transformation matrix in non-pivot rows
	mzd_t *R = mzd_init_window(A0, r-l, 0, A0->nrows, A0->ncols);
	mzd_set_ui(R, 0);

	mzd_free_window(A0);
	mzp_free(P);
	mzp_free(Q);
	mzd_free_window(A1);
	mzd_free_window(A00);
	mzd_free_window(A10);
	mzd_free_window(A01);
	mzd_free_window(A11);
	mzd_free_window(R);
	return 0;
}
*/


sw_list* isd(mzd_t* HzeroT, unsigned int l, unsigned int l2, unsigned int l3, unsigned int p, unsigned int e1, unsigned int e2, unsigned int w, unsigned int N, word** synds, unsigned int weight_threshold,unsigned int csize, unsigned long long max_iter, unsigned long long max_sol, unsigned long long max_time, ranctx* state, unsigned int skip) {
	unsigned int i, j ,ii;
	unsigned int n = HzeroT->nrows;
	unsigned int r = HzeroT->ncols;
	unsigned int k = n-r;

	unsigned int eff_word_len = min((unsigned int) HzeroT->ncols, (word_len*csize)); // number of bits to consider if we handle one word of data.

	mzd_t* HT  = mzd_init(   n,     r);
	mzd_t* BT  = mzd_init( k+l,     r);
	mzd_t* A   = mzd_init(   r,   r-l);
	mzd_t* AT  = mzd_init( r-l,     r);
	mzd_t* A_I = mzd_init(   r, 2*r-l);
	mzd_t* U   = mzd_init(   r,     r);

	mzd_t* Uprimemod  = mzd_init(eff_word_len, r);
	mzd_t* UprimemodT = mzd_init(r, eff_word_len);
	mzd_t* Usecondmod = mzd_init(r-eff_word_len, r);

	mzd_t* HprimemodT = mzd_init(k+l, eff_word_len);

	word* simple_HprimemodT = (word*) malloc((k+l)*sizeof(word)*csize); // copy of the internal data of HprimemodT; removes one layer of pointers and realigns datas

	unsigned int* perm = (unsigned int*) malloc(n*sizeof(unsigned int));
	unsigned int* perm_inv = (unsigned int*) malloc(n*sizeof(unsigned int));

	word* syndsprime;

	// If no syndrom is given, we assume the null vector as the target 
	// TODO : this will leak
	if (N == 0) {
		synds = (word**) malloc(sizeof(word*));
		synds[0] = (word*) calloc(bit_in_words(r), sizeof(word));
		N = 1;
	}
	
	/*
	 *  Syndsprime construction:
	 *  let X, Y, Z be different syndroms (N=3)
	 *  =>  Syndsprime = [X1 Y1 Z1 X2 Y2 Z2 ... Xcsize Ycsize Zcsize]
	 *  with X = [X1 X2 X3 ... Xcsize] each Xi of size 64 (or less if csize*64 > r)
	 */
	syndsprime = (word*) malloc(N*sizeof(word)*csize);

	sw_list* h = NULL;

	sub_isd_init(simple_HprimemodT, N, syndsprime, n, r, l, l2, l3, p, e1, e2, w, weight_threshold,csize, &h);
	final_test_init(r, w, BT, Usecondmod, synds);

	printf("n : %d\n", n);
	printf("r : %d\n", r);
	printf("w : %d\n", w);
	printf("N : %d\n", N);
	printf("l : %d\n", l);
	printf("l2 : %d\n", l2);
	printf("l3 : %d\n", l3);
	printf("p : %d\n", p);
	printf("e1 : %d\n", e1);
	printf("e2 : %d\n", e2);
	printf("max_iter : %lld\n", max_iter);
	printf("max_time : %lld\n", max_time);
	printf("max_sol : %lld\n", max_sol);
	printf("eff_word_len : %d\n", eff_word_len);
	printf("threshold : %d\n", weight_threshold);
	printf(" csize : %d\n", csize);

	for (i = 0; i < skip; ++i) {
		generate_permutation(perm, perm_inv, n, state);
	}

	/* Maybe unnecessary but not harmful; apply a reversible linear
	transformation to the columns of HzeroT to try to prevent eventual
	degenerate cases. */
	scramble(HzeroT, N, synds);

	unsigned long long nb_iter = 0;
	unsigned long long nb_sol = 0;
	unsigned long long start_date = time(NULL);
	unsigned long long start_cycles = cpucycles();

	// Main loop. Stops when 
	//    max_iter iterations are done 
	// or max_sol solutions have been found 
	// or max_time seconds elapsed.
	// Setting one of the three variable to 0 disable the corresponding stop condition.
	while (stop == 0) {
		if (print_status) {
			printf("nb_iter : %lld\n", nb_iter);
			print_status = 0;
		}

		pivot_probe_start();
		generate_permutation(perm, perm_inv, n, state);

		/*
		void xor_lines(word* res, int width, word* line1, word* line2){
			unsigned int i;
			for (i = 0; i < bit_in_words(width); ++i) {
				res[i]=line1[i] ^ line2[i];
			}
		}

		fprintf(stderr, "WARNING : USING TESTING SYNDROM\n");
		for (i = 0; i < bit_in_words(r); ++i) {
			synds[0][i] = 0;
		}

		// dumer 4 debug
		for (i = 0; i < w-4; ++i) {
			xor_lines(synds[0], r, synds[0], HzeroT->rows[i]);
		}
		for (i = 0; i < 2; ++i) {
			printf("%d : %lx\n", (r-l)+i, HzeroT->rows[(r-l)+i][0]);
			xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+i]);
		}
		for (i = 0; i < 2; ++i) {
			printf("%d : %lx\n", (r-l)+(k+l)/2+i, HzeroT->rows[(r-l)+(k+l)/2+i][0]);
			xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+(k+l)/2+i]);
		}
		printf("synd : %lx\n", synds[0][0]);

		// BJMM 8 debug

		for (i = 0; i < w-8; ++i) {
			xor_lines(synds[0], r, synds[0], HzeroT->rows[i]);
		}

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+(k+l)/2][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+(k+l)/2]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+4][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+4]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+5][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+5]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+1][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+1]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+(k+l)/4][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+(k+l)/4]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+2][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+2]);

		printf("%d : %lx\n", (r-l), HzeroT->rows[(r-l)+(k+l)/4+1][0]);
		xor_lines(synds[0], r, synds[0], HzeroT->rows[(r-l)+(k+l)/4+1]);

		printf("synd : %lx\n", synds[0][0]);

		generate_id_permutation(perm, perm_inv, n, state);

		*/


		apply_permutation(HT, HzeroT, perm, n);

		mzd_submatrix(BT, HT, r-l, 0, n, r);
		mzd_submatrix(AT, HT, 0, 0, r-l, r);
		mzd_transpose(A, AT);

		mzd_concat_id(A_I, A);

		// Abort this iteration if echelonization failed. We could have swapped
		// columns and adjusted the permutation but we prefer to abandon this
		// iteration and continue. 
		// This should happen with probability 2^-l
		
		if(mzd_partial_echelonize(A_I, l) == 1) {
			pivot_probe_stop();
			continue;
		}

		mzd_submatrix(U, A_I, 0, r-l, r, r+r-l);
		mzd_submatrix(Uprimemod, U, r-eff_word_len, 0, r, r);
		mzd_submatrix(Usecondmod, U, 0, 0, r-eff_word_len, r);
		mzd_transpose(UprimemodT, Uprimemod);
		mzd_mul_m4rm(HprimemodT, BT, UprimemodT, 0);

		//mzd_to_png(HprimemodT, "out.png", 0, NULL, 0);
		for (ii=0; ii<csize; ii++){
			for (i = 0; i < k+l; ++i) {
				simple_HprimemodT[i+ii*(k+l)] = HprimemodT->rows[i][ii];
			}
		}
		

		// Apply the transformation to the syndromes
		unsigned int s;
		unsigned int t;
		for (t=0; t< csize; t++){
			for (s = 0; s < N; ++s) {

				syndsprime[s+t*N] = 0;
				word chunk;

				for (i = 0; i < eff_word_len; ++i) {
					chunk = 0;
					for (j = 0; j < bit_in_words(r); ++j) {
						chunk ^= synds[s][j] & Uprimemod->rows[i][j];
					}
					if(isd_parity(chunk)){
						syndsprime[s+t*N] ^= (1UL << (i%64));
					}
				}
			}
		}

		pivot_probe_stop();
		bday_probe_start();
		sub_isd();
		bday_probe_stop();
		if (h) {
			sw_list* ptr = h;
			while(ptr) {
				++nb_sol;
				ptr = ptr->next;
			}
			process_solutions_on_the_fly(&h, w, l, BT, synds, U, perm_inv, nb_iter);
		}
		++nb_iter;

		if (max_iter != 0 && nb_iter >= max_iter) {
			stop = 1;
		}
		if (max_time != 0 && time(NULL) - start_date >= max_time) {
			stop = 1;
		}
		if (max_sol != 0 && nb_sol >= max_sol) {
			stop = 1;
		}
	}
	
	printf("Time spent : %llus\n", time(NULL) - start_date);
	long long pivot, bday, final;
	get_costs(nb_iter, &pivot, &bday, &final);
	sub_isd_report((cpucycles() - start_cycles)/nb_iter, pivot, bday, final);
	display(nb_iter);
	process_solutions_at_end(&h, w, l, BT, synds, U, perm_inv);

	sub_isd_free();
	final_test_free();

	mzd_free(HT);
	mzd_free(BT);
	mzd_free(A);
	mzd_free(AT);
	mzd_free(A_I);
	mzd_free(U);

	mzd_free(Uprimemod);
	mzd_free(UprimemodT);
	mzd_free(Usecondmod);

	mzd_free(HprimemodT);

	free(simple_HprimemodT);

	free(perm);
	free(perm_inv);
	free(syndsprime);
	return h;
}
