/**
 * \file process_solutions.h
 * \brief How to handle found solutions
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef PROCESS_SOLUTIONS_H
#define PROCESS_SOLUTIONS_H
#include <m4ri/m4ri.h>
#include "sparse_words_list.h"

/**
 * \brief Function called on each solution an iteration could find
 * \param partial_sol Sparse word of weight p giving a solution
 * \param w Weight of the search word
 * \param l Size of the window
 * \param BT Non permuted original matrix (r-l first rows truncated)
 * \param synds List of target syndromes
 * \param U Transition matrix that transformed HzeroT in partially eliminated form
 * \param perm_inv Inverse permutation of the current permutation
 * \param nb_iter Number of the current iteration
 * \todo Use isd_param struct
 */
void process_solution(sw_list* partial_sol, unsigned int w, unsigned int l, mzd_t* BT, word** synds, mzd_t* U, unsigned int* perm_inv, unsigned long long nb_iter);

#endif
