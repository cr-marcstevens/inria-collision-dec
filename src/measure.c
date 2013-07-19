#include <stdio.h>
#include "measure.h"
#include "cpucycles/cpucycles.h"

static long long pivot_cycles = 0;
static long long bday_cycles = 0;
static long long final_test_cycles = 0;
static long long nb_final_test = 0;
static long long nb_collision = 0;

void pivot_probe_start() {
	pivot_cycles -= cpucycles();
}
void pivot_probe_stop() {
	pivot_cycles += cpucycles();
}
void bday_probe_start() {
	bday_cycles -= cpucycles();
}
void bday_probe_stop() {
	bday_cycles += cpucycles();
}
void final_test_probe_start() {
	final_test_cycles -= cpucycles();
}
void final_test_probe_stop() {
	final_test_cycles += cpucycles();
}
void incr_final_test_counter() {
	++nb_final_test;
}
void incr_nb_collision_counter() {
	++nb_collision;
}

void get_costs(unsigned long long nb_iter, long long* pivot_cost, long long* bday_cost, long long* final_test_cost) {
	*pivot_cost = pivot_cycles/nb_iter;
	*bday_cost = bday_cycles/nb_iter;
	*final_test_cost = final_test_cycles/nb_final_test;
}

void display(unsigned long long nb_iter) {
		printf("\n");
		printf("total\n");
		
		printf("nb_iter : %lld\n", nb_iter);
		
		printf("collisions          %12lld collisions  \n", nb_collision);
		printf("test final count    %12lld tests       (%5.2f%%)\n", nb_final_test, 100.0*nb_final_test/nb_collision);
		long sum_cycles = pivot_cycles+bday_cycles;
		printf("pivot cycles        %12lld cycles      (%5.2f%%)\n", pivot_cycles, 100.0*pivot_cycles/sum_cycles);
		printf("bday cycles         %12lld cycles      (%5.2f%%)\n", bday_cycles - final_test_cycles, 100.0*(bday_cycles - final_test_cycles)/sum_cycles);
		if (nb_final_test != 0) {
			printf("test final cycles   %12lld cycles      (%5.2f%%)\n", final_test_cycles, 100.0*final_test_cycles/sum_cycles);
		}
		printf("\n");
		printf("avg on 1 iter\n");
		printf("pivot cycles        %12lld cycles/iter\n", pivot_cycles/nb_iter);
		printf("bday cycles         %12lld cycles/iter\n", bday_cycles/nb_iter);
		printf("collisions          %12lld collisions/iter\n", nb_collision/nb_iter);
		printf("\n");
}
