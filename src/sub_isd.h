#ifndef SUB_ISD_H
#define SUB_ISD_H
#include "isd.h"
#include "sparse_words_list.h"

void sub_isd_init(isd_params* params, word* local_L, word* local_synds, unsigned int local_N, sw_list** local_h, ranctx* state);
void sub_isd();
void sub_isd_report(unsigned long long cycles_per_iter, long long pivot_cost, long long bday_cost, long long final_test_cost);
void sub_isd_free();

#endif
