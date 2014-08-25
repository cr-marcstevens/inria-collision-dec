/**
 * \file isd.h
 * \brief Information Set Decoding
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef ISD_H
#define ISD_H
#include <m4ri/m4ri.h>
#include "sparse_words_list.h"
#include "prng.h"

/* A column index. This bounds k+l */
typedef uint16_t ci_t;

/**
 * \brief Parameters of the ISD algorithm and the sub_isd algorithm
 * \note Some parameters are specific to a sub_isd algorithm and won't impact others
 */
typedef struct {
	unsigned int n;                 /**< Length of the code (i.e. number of rows of HzeroT */
	unsigned int r;	                /**< Codimension of the code (i.e. number of columns of HzeroT */
	unsigned int k;	                /**< Dimension of the code (i.e. n-r)*/
	unsigned int l;                 /**< Optimisation parameter (all) */
	unsigned int l2;                /**< Optimisation parameter (MMT, BJMM) */
	unsigned int l3;                /**< Optimisation parameter (BJMM) */
	unsigned int alpha;             /**< Optimisation parameter (MMT) */
	unsigned int p;                 /**< Optimisation parameter (all) */
	unsigned int e1;                /**< Optimisation parameter (BJMM) */
	unsigned int e2;                /**< Optimisation parameter (BJMM) */
	unsigned int w;                 /**< Weight of the searched word */
	unsigned int weight_threshold;  /**< Optimisation parameter (all) see compute_threshold.c */
	unsigned int csize;             /**< BJMM specific parameter; number of word used to represent a column */

	unsigned long long max_iter;    /**< Stop condition : Stop when max_iter iterations are done */
	unsigned long long max_sol;     /**< Stop condition : Stop when max_sol solutions are found */
	unsigned long long max_time;    /**< Stop condition : Stop when max_time seconds have passed since launch */
} isd_params;

/**
 * \brief Print status. Triggered when receiving SIGUSR1
 */
void status_handler();

/**
 * \brief Finish current iteration and stop. Triggered when receiving SIGTERM or SIGINT
 */
void stop_handler();

/**
 * \brief Perform Information Set Decoding
 *
 * Main function of the program.
 * Search for params->w rows of HzeroT summing to one of the syndromes in synds
 * \param	HzeroT
 * \param N Number of syndromes in synds
 * \param synds Target syndromes
 * \param params isd and sub_isd algorithms parameters
 * \param state PRNG state
 * \param skip Number of iteration to skip before starting the search
 *
 * This algorithm :
 * 	- Permutes randomly columns of Hzero
 * 	- Perform a partial gaussian elimination on the permuted matrix (r-l columns processed)
 * 	- Extract the submatrix Hprime
 * 	- Call sub_isd on Hprime
 * 	- Process solutions output of sub_isd
 *
 * 	The function stops when :
 * 	 - params->max_iter iterations are done
 * 	 - params->max_sol solutions are found
 * 	 - params->max_time seconds have passed since launch
 * 	 - SIGTERM or SIGINT signal is received
 */
sw_list* isd(mzd_t* HzeroT, unsigned int N, word** synds, isd_params* params, ranctx* state, unsigned int skip);

#endif
