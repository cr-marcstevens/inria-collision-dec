/**
 * \file libisd.h
 * \brief Misc tools
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef LIBISD_H
#define LIBISD_H
#include <stdint.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>
#include "prng.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define word_size (sizeof(word))
#define word_len (word_size * 8)
#define bit_in_words(i) ((unsigned int) (((i) + word_len - 1) / word_len))
#define bit_in_bytes(i) ((unsigned int) (((i) + 8 - 1) / 8))

#define vec_set_bit(vec, j) ((vec)[(j)/word_len] |= (1L << ((j)%word_len)))
#define vec_get_bit(vec, j) (((vec)[(j)/word_len] >> ((j)%word_len)) &1)

/**
 * \brief Compute binomial coefficient n choose r
 */
double nCr(int n, int r);

/**
 * \brief Print binary representation of word (MSB to the left)
 */
void print_bin(word x);

/**
 * \brief Compare two words (usable with qsort)
 */
int word_cmp (const void * a, const void * b);

/**
 * \brief Compute scalar product a*b^t
 */
int pscal(word* a, word* b, int nb_bits);

/**
 * \brief Compute weight of a word (number of non zero position)
 */
unsigned int isd_weight(word w);

/**
 * \brief Compute parity of a word (1 if weight is odd; 0 otherwise)
 */
unsigned int isd_parity(word w);

#endif
