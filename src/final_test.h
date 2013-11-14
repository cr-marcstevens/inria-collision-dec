/**
 * \file final_test.h
 * \brief Final test of a candidate.
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#include <m4ri/m4ri.h>

/**
 * \brief Initialize the final_test() routine. Allocate required memory, put parameters and pointers to inputs in static global variables for future call to final_test().
 */
void final_test_init(unsigned int local_r, unsigned int local_w, mzd_t* local_BT, mzd_t* local_Usecondmod, word** local_synds);

/**
 * \brief Check whether a given p-uplet of columns indices output of sub_isd() gives a sum with weight lower than w-p on the complete matrix.
 * \param synd_idx The index of the syndrome considered
 * \param weight The weight of the sum of the p columns and considered syndrome on the lower part of the matrix
 * \param p The number of columns
 * \param ... The p columns
 */
int final_test(unsigned int synd_idx, unsigned int weight, unsigned int p, ... );

/**
 * \brief Like final_test but the p columns are passed as an array.
 */
int final_test_array(unsigned int synd_idx, unsigned int weight, unsigned int p, unsigned short* columns);

/**
 * \brief Free what final_test_init() allocated.
 */
void final_test_free();
