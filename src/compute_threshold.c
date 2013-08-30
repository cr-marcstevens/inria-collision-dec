/**
 * \file compute_threshold.c
 * \brief Compute the optimal threshold for the weight of the sum of p columns using exceeding bits, given relative cost of the different steps of the algorithm.
 * 
 * This program reads its parameters on stdin. See the sscanf() function call below to see the expected format. 
 * It computes the probability to trigger a final test given a random input sum (p_false_alarm)
 * and the probability not to trigger a final test given an input sum giving the solution (1 - p_not_miss).
 * Finally it minimizes (K0 + K1 + (p_false_alarm * K2 )) / p_not_miss where for one iteration K0 = cost of gaussian elimination; K1 = cost of subisd; and K2 = cost of the final tests
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#include <stdio.h>
#include <stdlib.h>
#include "libisd.h"

double combin_d(int n, int k){
	double x, y, i;
	if(k > n){
		return 0;
	}
	if (k > n/2){
		k = n-k;
	}
	x = 1;
	y = 1;
	i = n-k+1;
	while (i <= n){
		x = (x*i)/y;
		y += 1;
		i += 1;
	}
	return x;
}

int main()
{
	FILE* f = stdin;
	int i, T;
	int k=0, r=0, w=0, p=0, l=0, a=0;
	unsigned long long K0, K1, K2;
	double p_false_alarm = 0; /* Probability to trigger a final test given a random input sum */
	double p_not_miss = 0; /* Probability *not* to trigger a final test given an input sum giving the solution */
	char input[512];
	double cost;
	unsigned long long min_cost = 0;
	int arg_min = -1;
	int sucess = 0;
	do {
		if(fgets(input, 512, f) == NULL) {
			break;
		}
	} while ((sucess = sscanf(input, "couts pour k=%d r=%d w=%d p=%d l=%d : %lld %lld %lld", &k, &r, &w, &p, &l, &K0, &K1, &K2)) != 8 && (!feof(f)));


	int eff_word_bits = min(r, (int) word_len);
	a = eff_word_bits - l;

	if (sucess == 0) {
		fprintf(stderr, "error : no data found\n");
		exit(EXIT_FAILURE);
	}

	printf("k=%d\nr=%d\nw=%d\np=%d\nl=%d\na=%d\nK0=%lld\nK1=%lld\nK2=%lld\n",k,r,w,p,l,a,K0,K1,K2);
	printf("%2s : %12s %13s %40s\n", "T", "p_not_miss", "p_false_alarm", "cost");
	for (T = 1; T < a; ++T) {
		p_false_alarm = 0;
		p_not_miss = 0;
		for (i = 0; i < T; ++i) {
			p_false_alarm += combin_d(a, i);
		}
		p_false_alarm /= (1ULL << a);

		for (i = 0; i < T; ++i) {
			p_not_miss += combin_d(a, i)*combin_d(r-l-a, w-p-i);
		}
		p_not_miss /= combin_d(r-l, w-p);

		cost = (K0 + K1 + (p_false_alarm * K2 )) / p_not_miss;
		printf("%2d : %12g %13g %40f\n", T, p_not_miss, p_false_alarm, cost);
		if ((arg_min == -1) || (cost < min_cost)) {
			arg_min = T;
			min_cost = cost;
		}
	}
	printf("\n");

	if (K2 == 0) {
		fprintf(stderr, "Final test cost was 0, this value could be wrong\n");
	}
	printf("%d\n", arg_min);
	return 0;
}
