#ifndef ISD_H
#define ISD_H
#include "m4ri/m4ri.h"
#include "sparse_words_list.h"
#include "prng.h"

void status_handler();
void stop_handler();
sw_list* isd(mzd_t* HzeroT, unsigned int l, unsigned int l2, unsigned int l3, unsigned int p, unsigned int e1, unsigned int e2, unsigned int w, unsigned int N, word** synds, unsigned int weight_threshold, unsigned long long max_iter, unsigned long long max_sol, unsigned long long max_time, ranctx* state, unsigned int skip);

#endif
