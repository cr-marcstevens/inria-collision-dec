/**
 * \file io.h
 * \brief ISD input/output routines
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef IO_H
#define IO_H
#include <stdio.h>
#include <m4ri/m4ri.h>

/**
 * \brief Read chars from fd until char '0' or char '1' or EOF is read.
 * \return '0', '1' or EOF accordingly.
 */
int read_bit(FILE* fd);

/**
 * \brief Read a CSD problem instance from fd. Fill HzeroT, w, N and synds with data read.
 * \param HzeroT Input matrix
 * \param w Weight of the searched word
 * \param N Number of target syndromes
 * \param synds Target syndromes
 * \param fd Opened file descriptor
 * \note Exit execution on invalid file
 */
void CSD_from_file(mzd_t** HzeroT, unsigned int* w, unsigned int* N, word*** synds, FILE* fd);

/**
 * \brief Free memory allocated by CSD_from_file().
 */
void CSD_free(mzd_t* HzeroT, unsigned int N, word** synds);
#endif
