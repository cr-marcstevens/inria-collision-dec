#ifndef PROCESS_SOLUTIONS_H
#define PROCESS_SOLUTIONS_H
#include "m4ri/m4ri.h"
#include "sparse_words_list.h"

void process_solutions_on_the_fly(sw_list** eprime, unsigned int w, unsigned int l, mzd_t* BT, word** synds, mzd_t* U, unsigned int* perm_inv, unsigned long long nb_iter);
void process_solutions_at_end(sw_list** h, unsigned int w, unsigned int l, mzd_t* BT, word** synds, mzd_t* U, unsigned int* perm_inv);
#endif
