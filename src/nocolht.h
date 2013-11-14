/**
 * \file nocolht.h
 * \brief Hash table storing words discarding previous insertions in case of collision
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef NOCOLHT_H
#define NOCOLHT_H
#include <m4ri/m4ri.h>
typedef word nocolht_elt;
typedef nocolht_elt* nocolht;

/**
 * \brief Initialize an hash table of size s
 */
nocolht nocolht_init(unsigned long long size);

void nocolht_store(nocolht L, word index, nocolht_elt value);

/**
 * \brief Get first word contained at index i in L
 * \return A word if index i contains something; NONE otherwise
 */
nocolht_elt nocolht_get(nocolht L, word index);

/**
 * \brief Empty the table L.
 */
void nocolht_reset(nocolht L, unsigned long long L_size);

/**
 * \brief Free memory allocated by nocolht_init()
 */
void nocolht_free(nocolht L);
#endif
