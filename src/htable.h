#ifndef HTABLE_H
#define HTABLE_H

#define NONE -1

void htable_init(unsigned long long max_idx, unsigned long long nb_elt);
void htable_store(unsigned int index, int value);
int htable_get(unsigned int index);
int htable_next(unsigned int index, int value);
void htable_reset();
void htable_free();
void htable_stats();
#endif
