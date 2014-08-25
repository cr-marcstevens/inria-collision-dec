/**
 * \file support.h
 * \brief Functions used to cut a list in halves in differents ways.
 * This allow to restrict a list of columns to a support.
 *
 * Let L be the input list and Lr the output i.e. L restricted to the wanted support
 * Each "way" is composed of three functions :
 * 	prepare_* computes the length of Lr and allocate the coresponding memory
 * 	build_* effectively copy datas from L to Lr
 * 	inv_* gives back the index in L of the index of an element in Lr
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef SUPPORT_H
#define SUPPORT_H
#include <m4ri/m4ri.h>
#include "isd.h"

/* First half */
void prepare_half0(word** Lr, unsigned int* Lr_len, unsigned int L_len);
void build_half0(word* Lr, word* L, unsigned int L_len);
ci_t inv_half0(ci_t i, unsigned int L_len);

/* Second half */
void prepare_half1(word** Lr, unsigned int* Lr_len, unsigned int L_len);
void build_half1(word* Lr, word* L, unsigned int L_len);
ci_t inv_half1(ci_t i, unsigned int L_len);

/* First and third quarter */
void prepare_quarter0and2(word** Lr, unsigned int* Lr_len, unsigned int L_len);
void build_quarter0and2(word* Lr, word* L, unsigned int L_len);
ci_t inv_quarter0and2(ci_t i, unsigned int L_len);

/* Second and fourth quarter */
void prepare_quarter1and3(word** Lr, unsigned int* Lr_len, unsigned int L_len);
void build_quarter1and3(word* Lr, word* L, unsigned int L_len);
ci_t inv_quarter1and3(ci_t i, unsigned int L_len);
#endif
