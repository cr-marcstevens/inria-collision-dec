#ifndef IAHT_H
#define IAHT_H
#include <stdint.h>
typedef uint32_t iaht_elt;
typedef iaht_elt* iaht_list;
typedef iaht_list* iaht;

iaht iaht_init(unsigned long long size);
void iaht_store(iaht L, word index, int x, ...);
iaht_list iaht_get(iaht L, word index);
void iaht_reset(iaht L, unsigned long long size);
void iaht_free(iaht L, int size);
#endif
