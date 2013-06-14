#ifndef SUB_ISD_H
#define SUB_ISD_H
#include "libisd.h"
#include "sparse_words_list.h"

extern unsigned int p;

void sub_isd_init(word* simple_HprimemodT, unsigned int N, word* syndprime, unsigned int local_n, unsigned int local_r, unsigned int local_l, unsigned int local_l2, unsigned int local_w, unsigned int local_threshold, sw_list** h);
void sub_isd();
void sub_isd_report(unsigned long long cycles_per_iter);
void sub_isd_free();

#endif
