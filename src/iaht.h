/**
 * \file iaht.h
 * \brief Array hash table containing list of integer x-uplet
 *
 * First element of each lists is the number of x-uplet in the list. 
 * Next x elements are the first x-uplet and so on...
 * Lists are reallocated when an element is added.
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */

#ifndef IAHT_H
#define IAHT_H
#include <stdint.h>
typedef uint32_t iaht_elt;
typedef iaht_elt* iaht_list;
typedef iaht_list* iaht;

/**
 * \brief Initialize a hash table of size s
 */
iaht iaht_init(unsigned long long s);

/**
 * \brief Store x-uplet at index i of L.
 * \param L The hash table
 * \param i The index
 * \param x number of element in the x-uplet
 * \param ... x elements
 */
void iaht_store(iaht L, word i, int x, ...);

/**
 * \brief Return pointer to the first list of x-uplet at index i
 */
iaht_list iaht_get(iaht L, word i);

/**
 * \brief Return next list given that ptr was returned previously (by iaht_get() or iaht_next()) and p elements where used.
 */
iaht_list iaht_next(iaht L, word i, iaht_list ptr, unsigned int p);

/**
 * \brief Empty the table L
 */
void iaht_reset(iaht L, unsigned long long size);

/**
 * \brief Free memory allocated by iaht_init().
 */
void iaht_free(iaht L, unsigned long long size);
#endif
