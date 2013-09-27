/**
 * \file counterht_col2.h
 * \brief Hash table storing counters
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef COUNTERHT_H
#define COUNTERHT_H
#include <stdint.h>
#include "libisd.h"

#define NONE -1
typedef int32_t counter;
typedef counter counter_container;

typedef struct {
	counter_container* main_table;
	counter_container* col_table;
	unsigned long long main_table_size;
	unsigned long long col_table_size;
} *counterht;

/**
 * \brief Return the counter contained in the container
 */

counter counter_container_open(counter_container* c);

/**
 * \brief Initialize an hash table of size s
 */
counterht counterht_init(unsigned long long nb_of_hash_val, unsigned long long nb_of_insert);

void counterht_store(counterht L, word index, counter value);

/**
 * \brief Get first counter contained at index i in L
 * \return A counter if index i contains something; NONE otherwise
 */
counter_container* counterht_get(counterht L, word index);

/**
 * \brief Get next counter contained at index i in L given that we got counter c during a previous call to counterht_next or counterht_get
 */
counter_container* counterht_next(counterht L, word i, counter_container* c);

/**
 * \brief Empty L
 */
void counterht_reset(counterht L, unsigned long long L_size);
void counterht_free(counterht L);

#endif
