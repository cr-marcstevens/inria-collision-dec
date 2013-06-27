#include <stdlib.h>
#include "sub_isd.h"
#include "m4ri/m4ri.h"
#include "libisd.h"
#include "final_test.h"
#include "htable.h"
#include "sparse_words_list.h"
#include "time.h"
#include "measure.h"
#include "cpucycles/cpucycles.h"

#define p 6
static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, w, L_len, threshold;
static int shift;
static word* xors_table;
static sw_list** h;

void unpack3_counter(unsigned int* c1, unsigned int* c2, unsigned int* c3, unsigned int counter) {
	unsigned int col_1 = exp(log(6*counter)/3);
	do {
		++col_1;
	} while (counter > (col_1*(col_1-2)*(2*col_1-2)/12));
	--col_1;
	unsigned int col_2 = (int) (sqrt(2.0f*(counter-(col_1*(col_1-2)*(2*col_1-2)/12))+0.25f) + 0.5f);
	unsigned int col_3 = counter - (col_1*(col_1-2)*(2*col_1-2)/12) - (col_2*(col_2-1))/2;
	*c1 = col_1; *c2 = col_2; *c3 = col_3;
}

void sub_isd_init(word* simple_HprimemodT, unsigned int local_N, word* local_syndsprime, unsigned int local_n, unsigned int local_r,unsigned int local_l, unsigned int local_l2, unsigned int local_l3, unsigned int local_p, unsigned int local_e1, unsigned int local_e2, unsigned int local_w, unsigned int local_threshold,unsigned int local_csize, sw_list** local_h) {
	L = simple_HprimemodT;
	N = local_N;
	syndsprime = local_syndsprime;
	n = local_n;
	r = local_r;
	l = local_l;
	(void) local_l2;
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
	shift = min(r, word_len) - l;

	htable_init(1ULL << l, nCr(L_len, p/2));
	xors_table = (word*) malloc(nCr(L_len, p/2) * sizeof(word));
}

void sub_isd() {
	word synd = syndsprime[0]; //DOOM not implemented

	unsigned int c1, c2, c3, c4, c5, c6;
	unsigned int counter;
	int current;
	word res;
	unsigned int weight;
	int final_weight;

	htable_reset();
	counter = 0;
	for(c1 = 2; c1 < L_len/2; ++c1) {
		for(c2 = 1; c2 < c1; ++c2) {
			for(c3 = 0; c3 < c2; ++c3) {
				res = L[c1] ^ L[c2] ^ L[c3];
				htable_store(res >> shift, counter);
				xors_table[counter] = res;
				++counter;
			}
		}
	}
	for (c4 = L_len/2 + 2; c4 < L_len; ++c4) {
		for (c5 = L_len/2 + 1; c5 < c4; ++c5) {
			for (c6 = L_len/2; c6 < c5; ++c6) {
				res = synd ^ L[c4] ^ L[c5] ^ L[c6];
				for(current = htable_get(res >> shift); current != NONE; current = htable_next(res >> shift, current)) {
				incr_nb_collision_counter();
					weight = isd_weight(xors_table[current] ^ res);
					if (weight <= threshold) {
						incr_final_test_counter();
						final_test_probe_start();
						unpack3_counter(&c1, &c2, &c3, current);
						final_weight = final_test(0, weight, p, c1, c2, c3, c4, c5, c6);
						final_test_probe_stop();
						if (final_weight != -1) {
							*h = sw_list_add(*h, 0, final_weight, p, c1, c2, c3, c4, c5, c6);
						}
					}	
				}
			}
		}
	}
}

void sub_isd_report(unsigned long long cycles_periter) {

 /* p_miss is the proportion of candidates that are eliminated by the
 weight threshold condition whereas they where the solution. That's why
 we have to multiply the number of collision needed by this proportion
 */
	double p_miss = 0;
	unsigned int i;
	unsigned int eff_word_len = min(r, word_len);
	for (i = 0; i < threshold; ++i) {
		p_miss += nCr(eff_word_len - l, i)*nCr(r - eff_word_len, w-p-i);
	}
	p_miss /= nCr(r-l, w-p);
	p_miss = 1 - p_miss;
	printf("Proba to miss the solution : %8.3g %%\n", 100*p_miss);

	double nb_col_needed = nCr(n, w) / nCr(r-l, w-p) / (1ULL<<l);
	nb_col_needed += p_miss * nb_col_needed;
	double avg_nb_col_periter = nCr(L_len/2, p/2) * nCr(L_len - L_len/2, p/2) / (1ULL<<l);
	
	/*
	printf("nb_col_needed : %12.4g\n", (double)nb_col_needed);
	printf("avg_nb_col_periter : %12.4g\n", avg_nb_col_periter);
	double nb_iter_needed = nb_col_needed / avg_nb_col_periter;
	printf("nb_iter_needed : %12.4g\n", nb_iter_needed);
	printf("cycles_periter : %lld\n", cycles_periter);
	double cycles_needed = nb_iter_needed*cycles_periter;
	printf("cycles_needed (log2) : %12.4g\n", (double)log(cycles_needed)/log(2));
	printf("cpucycles_persecond() : %lld\n", cpucycles_persecond());
	*/
	double time_needed = (cycles_periter * nb_col_needed / avg_nb_col_periter) / cpucycles_persecond();
	
	printf("time needed (%.1fGHz) : %12.4gs : ", cpucycles_persecond()/1000000000.0, time_needed);
	time_t time = (time_t) time_needed;
	struct tm * tm_now = gmtime (&time);
	if (tm_now == NULL) {
		printf("more than 2^%ld years\n", 8*sizeof(int));
	}
	else {
		tm_now->tm_year -= 1970;
		tm_now->tm_mon -= 1;
		tm_now->tm_mday -= 1;
		//	tm_now->tm_hour -= 1;
		char s_now[80];
		strftime (s_now, sizeof s_now, "%Yy %mM %dd %Hh %Mm %Ss", tm_now);
		printf("%s\n", s_now);
	}
}

void sub_isd_free() {
	htable_free();
	free(xors_table);
}
