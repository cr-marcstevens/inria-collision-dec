#ifndef IO_H
#define IO_H
#include <stdio.h>
#include <m4ri/m4ri.h>
int read_bit(FILE* fd);
void CSD_from_file(mzd_t** HzeroT, unsigned int* w, unsigned int* N, word*** synds, FILE* fd);
void CSD_free(mzd_t* HzeroT, unsigned int N, word** synds);
#endif
