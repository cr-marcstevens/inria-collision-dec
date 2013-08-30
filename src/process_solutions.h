/**
 * \file process_solutions.h
 * \brief How to handle found solutions
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef PROCESS_SOLUTIONS_H
#define PROCESS_SOLUTIONS_H
#include "m4ri/m4ri.h"
#include "sparse_words_list.h"

/**
 * \brief Function called after each iteration that found one or more solutions
 * \param eprime Linked list of sparse word of weight p giving a solution
 * \param w Weight of the search word
 * \param l Size of the window
 * \param BT Non permuted original matrix (r-l first rows truncated)
 * \param synds List of target syndromes
 * \param U Transition matrix that transformed HzeroT in partially eliminated form
 * \param perm_inv Inverse permutation of the current permutation
 * \param nb_iter Number of the current iteration
 * \todo Use isd_param struct
 * \todo change it so that it works only on one solution instead of a list (and change subisds so that they filter there list of solutions before returning)
 */
void process_solutions_on_the_fly(sw_list** eprime, unsigned int w, unsigned int l, mzd_t* BT, word** synds, mzd_t* U, unsigned int* perm_inv, unsigned long long nb_iter);

/**
 * \brief Function called when the main loop stopped
 * \copydetails process_solutions_on_the_fly()
 */
void process_solutions_at_end(sw_list** h, unsigned int w, unsigned int l, mzd_t* BT, word** synds, mzd_t* U, unsigned int* perm_inv, unsigned long long nb_iter);
#endif
