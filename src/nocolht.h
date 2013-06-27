#ifndef NOCOLHT_H
#define NOCOLHT_H
#include "m4ri/m4ri.h"
typedef word nocolht_elt;
typedef nocolht_elt* nocolht;

nocolht nocolht_init(unsigned long long size);
void nocolht_store(nocolht L, word index, nocolht_elt value);
nocolht_elt nocolht_get(nocolht L, word index);
void nocolht_reset(nocolht L, unsigned long long L_size);
void nocolht_free(nocolht L);
#endif
