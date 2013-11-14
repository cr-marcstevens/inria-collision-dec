/**
 * \file waht.h
 * \brief Array hash table containing lists of words.
 *
 * First element of each lists contains the number of element of the list.
 * Lists are reallocated when an element is added.
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */

#ifndef WAHT_H
#define WAHT_H
#include <m4ri/m4ri.h>
typedef word* waht_list;
typedef waht_list* waht;

/**
 * \brief Initialize a hash table of size s
 */
waht waht_init(unsigned long long s);

/**
 * \brief Store word w at index i of L.
 */
void waht_store(waht L, word i, word w);

/**
 * \brief Return pointer to the first word stored at index i of L;
 * \return NULL if index i of L is empty; pointer to a word otherwise
 */
word* waht_get(waht L, word i);

/**
 * \brief Return next pointer to the word stored at index i of L given that a previous call to waht_get() or waht_next() returned prev.
 * \return NULL if index i of L is empty; pointer to a word otherwise
 */
word* waht_next(waht L, word index, waht_list prev);

/**
 * \brief Empty the table L
 */
void waht_reset(waht L, unsigned long long size);

/**
 * \brief Free memory allocated by waht_init()
 */
void waht_free(waht L, int size);
#endif
