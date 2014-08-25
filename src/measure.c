#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include "measure.h"
#include "cpucycles/cpucycles.h"
#include "isd.h"
#include "libisd.h"

static long long total_time = 0;

static long long total_cycles = 0;
static long long pivot_cycles = 0;
static long long bday_cycles = 0;
static long long final_test_cycles = 0;

static long long nb_iter = 0;
static long long nb_final_test = 0;
static long long nb_collision = 0;

void total_time_stopwatch_start() {
	total_time -= time(NULL);
}
void total_time_stopwatch_stop() {
	total_time += time(NULL);
}

void total_cycle_stopwatch_start() {
	total_cycles -= cpucycles();
}
void total_cycle_stopwatch_stop() {
	total_cycles += cpucycles();
}
void pivot_cycle_stopwatch_start() {
	pivot_cycles -= cpucycles();
}
void pivot_cycle_stopwatch_stop() {
	pivot_cycles += cpucycles();
}
void bday_cycle_stopwatch_start() {
	bday_cycles -= cpucycles();
}
void bday_cycle_stopwatch_stop() {
	bday_cycles += cpucycles();
}
void final_test_cycle_stopwatch_start() {
	final_test_cycles -= cpucycles();
}
void final_test_cycle_stopwatch_stop() {
	final_test_cycles += cpucycles();
}

void incr_iter_counter() {
	++nb_iter;
}
void incr_final_test_counter() {
	++nb_final_test;
}
void incr_collision_counter() {
	++nb_collision;
}

void reset() {
	total_time = 0;
	total_cycles = 0;
	pivot_cycles = 0;
	bday_cycles = 0;
	final_test_cycles = 0;
	nb_iter = 0;
	nb_final_test = 0;
	nb_collision = 0;
}

void report(isd_params* params) {
	if (nb_iter == 0) {
		fprintf(stderr, "No iteration done\n");
	}
	else {
		unsigned long long pivot_cost = pivot_cycles/nb_iter;
		unsigned long long bday_cost = bday_cycles/nb_iter;
		unsigned long long final_test_cost = (nb_final_test == 0) ? 0 : final_test_cycles/nb_final_test;
		unsigned long long cycles_periter = total_cycles / nb_iter;
		float iter_persecond;

		double p_miss = 0;
		unsigned int i;
		unsigned int d = min(params->r, word_len);

		double nb_col_needed;
		double nb_col_periter;
		double nb_iter_needed;
		double cycles_needed;
		double time_needed;
		time_t time;
		struct tm* tm_now;
		char s_now[80];

		unsigned int n, r, k, l, w, weight_threshold, p;

		n = params->n;
		r = params->r;
		k = params->k;
		l = params->l;
		w = params->w;
		p = params->p;
		weight_threshold = params->weight_threshold;

		/* If stdout is not a tty, we assume it is piped to compute_threshold; you can pipe to 'cat' to see the line */
		if (!isatty(fileno(stdout))) {
			printf("couts pour k=%d r=%d w=%d p=%d l=%d : %lld %lld %lld\n", params->k, params->r, params->w, params->p, params->l, pivot_cost, bday_cost, final_test_cost); // this line can be parsed by compute_threshold
		}

		if (total_time != 0) {
			iter_persecond = (float)nb_iter/total_time;
		}
		else {
			iter_persecond = (float)cpucycles_persecond() / cycles_periter;
		}

		printf("%lld iterations done in %lld seconds (%.2f iter/s)\n", nb_iter, total_time, iter_persecond);

		if (nb_iter == 0) {
			printf("No iteration done\n");
			return;
		}
		printf("\n");
		printf("Total cycles : %lld\n\n", total_cycles);
		printf("Per iteration : \n");

		printf("\t(k+l choose p)/2^l %13.2f\n", (float)nCr(k+l, p)/(1ULL<<l));
		printf("\tCollisions         %13.2f (%5.2f%%)\n", (float)nb_collision/nb_iter, 100*((float)nb_collision/nb_iter)/(nCr(k+l, p)/(1ULL<<l)));
		printf("\tFinal tests        %13.2f", (float)nb_final_test/nb_iter);
		if (nb_collision != 0) {
			printf(" (%5.2f%%)", 100*(float)nb_final_test/nb_collision);
		}
		printf("\n");
		printf("\n");

		printf("\tPivot cycles       %13lld (%5.2f%%)\n", pivot_cycles/nb_iter, 100.0*pivot_cycles/total_cycles);
		printf("\tBirthday cycles    %13lld (%5.2f%%)\n", bday_cycles/nb_iter, 100.0*bday_cycles/total_cycles);
		printf("\tFinal test cycles  %13lld", final_test_cycles/nb_iter);
		if (nb_collision != 0) {
			printf(" (%5.2f%%)", 100*(float)final_test_cycles/total_cycles);
		}
		printf("\n");
		printf("\n");


		/* p_miss is the proportion of candidates that are eliminated by the
			 weight threshold condition whereas they where the solution. That's why
			 we have to multiply the number of collision needed by this proportion
			 */

		for (i = 0; i < weight_threshold; ++i) {
			p_miss += nCr(d - l, i)*nCr(r - d, w-p-i);
		}
		p_miss /= nCr(r-l, w-p);

		nb_col_needed = nCr(n, w) / nCr(r-l, w-p) / (1ULL<<l);
		nb_col_needed += nb_col_needed / p_miss;
		//nb_col_periter = nCr((k+l)/2, p/2) * nCr((k+l) - (k+l)/2, p/2) / (1ULL<<l); // <- theoretical for dumer disjoint support
		nb_col_periter = nb_collision/nb_iter;
		nb_iter_needed = nb_col_needed / nb_col_periter;
		cycles_needed = nb_iter_needed*cycles_periter;
		time_needed = (cycles_periter * nb_col_needed / nb_col_periter) / cpucycles_persecond();

		time = (time_t) time_needed;
		tm_now = gmtime (&time);

		printf("Threshold : %u\n", params->weight_threshold);
		printf("Miss prob : %g\n", 1-p_miss);

		printf("Average requirement per solution\n");
		printf("\tCollisions (log2) : %12.4g\n", (double)log(nb_col_needed)/log(2));
		printf("\tIterations :        %12.4g\n", nb_iter_needed);
		printf("\tCycles (log2) :     %12.4g\n", (double)log(cycles_needed)/log(2));

		printf("\tTime (%.1fGHz) :    %12.4gs", cpucycles_persecond()/1000000000.0, time_needed);
		if (tm_now == NULL) {
			printf(" (more than 2^%ld years)", 8*sizeof(int));
		}
		else {
			tm_now->tm_year -= 1970;
			tm_now->tm_mon -= 1;
			tm_now->tm_mday -= 1;
			strftime (s_now, sizeof s_now, "%Yy %mM %dd %Hh %Mm %Ss", tm_now);
			printf("(%s)", s_now);
		}
		printf("\n");
	}
}
