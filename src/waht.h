#ifndef WAHT_H
#define WAHT_H
#include "m4ri/m4ri.h"
typedef word waht_elt;
typedef waht_elt* waht_list;
typedef waht_list* waht;

waht waht_init(unsigned long long size);
void waht_store(waht L, word index, waht_elt value);
waht_list waht_get(waht L, word index);
void waht_reset(waht L, unsigned long long size);
void waht_free(waht L, int size);
#endif
