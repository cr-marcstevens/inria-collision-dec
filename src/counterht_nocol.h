/**
 * \file counterht_nocol.h
 * \brief Hash table storing counters discarding previous insertions in case of collision
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef COUNTERHT_H
#define COUNTERHT_H
#include <stdint.h>
#include "libisd.h"

#define NONE -1
typedef int32_t counter;
typedef counter* counterht; 

/**
 * \brief Initialize an hash table of size s
 */
counterht counterht_init(unsigned long long s);

void counterht_store(counterht L, word index, counter value);

/**
 * \brief Get first counter contained at index i in L
 * \return A counter if index i contains something; NONE otherwise
 */
counter counterht_get(counterht L, word i);

/**
 * \brief Get next counter contained at index i in L given that we got counter c during a previous call to counterht_next or counterht_get
 * \note Since this hash table does not consider duplicate, this function always returns NONE; this function is here for future ways of storing counters
 */
counter counterht_next(counterht L, word i, counter c);

/**
 * \brief Empty L
 */
void counterht_reset(counterht L, unsigned long long L_size);
void counterht_free(counterht L);

#endif
