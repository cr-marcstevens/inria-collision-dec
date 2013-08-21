#include <stdlib.h>
#include <time.h>
#include "sub_isd.h"
#include "m4ri/m4ri.h"
#include "libisd.h"
#include "final_test.h"
#include "sparse_words_list.h"
#include "measure.h"
#include "cpucycles/cpucycles.h"

#define p 4
static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, w, L_len, threshold;
static unsigned int lprime;
static int shift;

#define NONE -1
static unsigned int L0_size;
static int* L0;
static word* xors_table;

static sw_list** h;

void unpack2_counter(unsigned int* c1, unsigned int* c2, unsigned int counter) {
	(*c1) = (int) (sqrt(2.0f*counter+0.25f) + 0.5f);
	(*c2) = counter - (*c1*(*c1-1))/2;
}

void print_parameters(isd_params* params) {
	printf("n : %d\n", params->n);
	printf("r : %d\n", params->r);
	printf("w : %d\n", params->w);
	printf("l : %d\n", params->l);
	printf("p : %d\n", p);
	printf("eff_word_len : %ld\n", min(params->r, word_len));
	printf("threshold : %d\n", params->weight_threshold);
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
	w = params->w;

	k = n-r;

	L_len = k+l;

	threshold = params->weight_threshold;

	unsigned long long nb_of_sums = nCr(L_len/2, p/2);
	lprime = (l+ log(nb_of_sums)/log(2))/2;
	printf("lprime : %d\n", lprime);

	shift = min(r, word_len) - lprime;

	L0_size = 1ULL << lprime;
	L0 = (int*) malloc(L0_size*sizeof(int));

	xors_table = (word*) malloc(nb_of_sums * sizeof(word));
}


void sub_isd() {
	word synd = syndsprime[0]; //DOOM not implemented

	unsigned int c1, c2, c3, c4;
	unsigned int counter;
	int current;
	word res;
	word sumS3; /* intermediate sums */
	unsigned int weight;
	int final_weight;
	word max_word_zero_l_bits = 1UL << (word_len - l); /* A word with its l MSB zeroed is lower than this value */

	memset(L0, NONE, L0_size*sizeof(*L0));
	counter = 0;
	for(c1 = 1; c1 < L_len/2; ++c1) {
		for(c2 = 0; c2 < c1; ++c2) {
			res = L[c1] ^ L[c2];
			L0[res >> shift] = counter;
			xors_table[counter] = res;
			++counter;
		}
	}
	/*
		 htable_stats();
		 exit(1);
		 */
	for (c3 = L_len/2 + 1; c3 < L_len; ++c3) {
		sumS3 = synd ^ L[c3];
		for (c4 = L_len/2; c4 < c3; ++c4) {
			res = sumS3 ^ L[c4];
			current = L0[res >> shift];
			if (current != NONE) {
				res ^= xors_table[current];
				if (res < max_word_zero_l_bits) {
					incr_nb_collision_counter();
					weight = isd_weight(res);
					if (weight <= threshold) {
						incr_final_test_counter();
						final_test_probe_start();
						unpack2_counter(&c1, &c2, current);
						final_weight = final_test(0, weight, p, c1, c2, c3, c4);
						final_test_probe_stop();
						if (final_weight != -1) {
							*h = sw_list_add(*h, 0, final_weight, p, c1, c2, c3, c4);
						}
					}
				}
			}
		}
	}
}

void sub_isd_report(unsigned long long cycles_periter, long long pivot_cost, long long bday_cost, long long final_test_cost) {

	printf("couts pour k=%d r=%d w=%d p=%d l=%d : %lld %lld %lld\n", k, r, w, p, l, pivot_cost, bday_cost, final_test_cost); // this line can be parsed by compute_threshold

	unsigned long long iter_persecond = cpucycles_persecond() / cycles_periter;
	printf("Iterations/secondes : %lld\n", iter_persecond);

 /* p_miss is the proportion of candidates that are eliminated by the
 weight threshold condition whereas they where the solution. That's why
 we have to multiply the number of collision needed by this proportion
 */
	unsigned int i;
	unsigned int eff_word_len = min(r, word_len);

	double nb_col_needed;
	double nb_col_periter;
	double nb_iter_needed;
	double cycles_needed;
	double time_needed;
	time_t time;
	struct tm* tm_now;
	char s_now[80];

	double p_miss = 0;
	for (i = 0; i < threshold; ++i) {
		p_miss += nCr(eff_word_len - l, i)*nCr(r - eff_word_len, w-p-i);
	}
	p_miss /= nCr(r-l, w-p);
	p_miss = 1 - p_miss;

	nb_col_needed = nCr(n, w) / nCr(r-l, w-p) / (1ULL<<l);
	nb_col_needed += p_miss * nb_col_needed;
	nb_col_periter = nCr(L_len/2, p/2) * nCr(L_len - L_len/2, p/2) / (1ULL<<l);
	nb_iter_needed = nb_col_needed / nb_col_periter;
	cycles_needed = nb_iter_needed*cycles_periter;
	time_needed = (cycles_periter * nb_col_needed / nb_col_periter) / cpucycles_persecond();

	time = (time_t) time_needed;
	tm_now = gmtime (&time);

	printf("nb_col_needed : %12.4g\n", (double)nb_col_needed);
	printf("nb_col_periter : %12.4g\n", nb_col_periter);
	printf("nb_iter_needed : %12.4g\n", nb_iter_needed);
	printf("cycles_periter : %lld\n", cycles_periter);
	printf("cycles_needed (log2) : %12.4g\n", (double)log(cycles_needed)/log(2));

	printf("collisions needed %12.4g\n", nb_col_needed);


	printf("time needed (%.1fGHz) : %12.4gs : ", cpucycles_persecond()/1000000000.0, time_needed);
	if (tm_now == NULL) {
		printf("more than 2^%ld years\n", 8*sizeof(int));
	}
	else {
		tm_now->tm_year -= 1970;
		tm_now->tm_mon -= 1;
		tm_now->tm_mday -= 1;
		strftime (s_now, sizeof s_now, "%Yy %mM %dd %Hh %Mm %Ss", tm_now);
		printf("%s\n", s_now);
	}
}

void sub_isd_free() {
	free(L0);
	free(xors_table);
}
