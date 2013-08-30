/**
 * \file sub_isd.h
 * \brief sub_isd algorithm
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */

#ifndef SUB_ISD_H
#define SUB_ISD_H
#include "isd.h"
#include "sparse_words_list.h"

/**
 * \brief Initialize the sub_isd() routine. Allocate required memory, put parameters and pointers to inputs in static global variables for future call to sub_isd().
 */
void sub_isd_init(isd_params* params, word* local_L, word* local_synds, unsigned int local_N, sw_list** local_h, ranctx* state);

/**
 * \brief Run one iteration of the sub_isd_routine. It gets its parameters and inputs from static global variables set by the sub_isd_init function.
 */
void sub_isd();

/**
 * \brief Free memory allocated by sub_isd_init().
 */
void sub_isd_free();

#endif
