/**
 * \file sparse_words_list.h
 * \brief Linked lists of sparse words
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef SPARSE_WORDS_LIST_H
#define SPARSE_WORDS_LIST_H

/** 
 * \brief Sparse word node
 */
typedef struct sw {
	unsigned int weight;      /**< Weight of the word */
	unsigned int *pos;        /**< Array of size weight containing the indices of the non-null positions of the word */
	unsigned int synd_idx;    /**< Index of the target syndrom considered */
	unsigned int synd_weight; /**< Weight of the syndrom of the word added to the target syndrom */
	unsigned int sorted;      /**< 1 if pos is sorted, 0 instead */

	struct sw *next;          /**< Next node in the linked list */
} sw;

typedef sw sw_list;

/**
 * \brief Allocate a new sparse word of weight p
 */
sw* sw_new(unsigned int p);
/**
 * \brief	Allocate a new sparse word of weight p and fill its fields with the input values
 * \param ... A va_list of p unsigned int with the non-null positions of the word
 */
sw* sw_filled_new(unsigned int synd_idx, unsigned int synd_weight, unsigned int p, ...);

/**
 * \brief Allocate a new sparse word of weight p and fill its fields with the input values
 * \param columns An array containing the p non-null positions of the word
 */
sw* sw_filled_new_array(unsigned int synd_idx, unsigned int synd_weight, unsigned int p, unsigned short* columns);
/**
 * \brief Sort the pos array
 */
void sw_sort(sw* h);
void sw_print(sw* h);
/**
 * \brief Compare two sparse words
 * \return 0 if w1 and w2 are identical; -1 otherwise
 * \note will call sw_sort() on both inputs
 */
int sw_cmp(sw_list* w1, sw_list* w2);
void sw_free(sw_list* h);

/**
 * \brief Allocate a new sparse word of weight p, fill its fields with the input values, fill the pos field with the value in the columns array and append it to h
 * TODO : delete this function and use sw_list_append(sw_list_new_array
 */
void sw_list_add_array(sw_list** h, unsigned int synd_idx, unsigned int synd_weight, unsigned int p, unsigned short* columns);

/**
 * \brief Append a sparse word to a sparse word list
 */
void sw_list_append(sw_list** h, sw* new);

/**
 * \brief Filter a sparse word list by removing duplicates
 * \note will call sw_sort() on each word in input list
 */
void sw_list_uniq(sw_list** h);
void sw_list_print(sw_list* h);

/**
 * \brief Free all words in a sparse word list
 */
void sw_list_free(sw_list* h);

#endif
