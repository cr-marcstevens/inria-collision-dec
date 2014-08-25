/**
 * \file ciht_t.h
 * \brief Hash table storing p-uplet of columns indices
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef CIHT_H
#define CIHT_H
#include <stdint.h>
#include "isd.h"
#include "libisd.h"

// A columns indices hash table. It is an array of p*nb_of_hash_val ci_t.
typedef ci_t* ciht_t;

// The value returned when no data is found at a given index in the hash table
#define NONE ((ci_t)-1)

/**
 * \brief Initialize an hash table of size nb_of_hash_val, that will store a maximum of nb_of_insert p-uple of column indices.
 */
ciht_t ciht_init(unsigned long long nb_of_hash_val, unsigned int p, unsigned long long nb_of_insert);

/**
	* \brief Store the p contiguous columns indices in L at index i
	*/
void ciht_store(ciht_t L, word i, ci_t* value, unsigned int p);

/**
 * \brief Get p contiguous ci_t stored at index i in L.
 * \return A pointer to p contiguous ci_t. These will all be NONE if nothing is stored at this index.
 */
ci_t* ciht_get(ciht_t L, word index, unsigned int p);

/**
 * \brief Empty L
 */
void ciht_reset(ciht_t L, unsigned long long L_size, unsigned int p);

/**
 * \brief Free memory. Undo what ciht_init did.
 */
void ciht_free(ciht_t L);

#endif
