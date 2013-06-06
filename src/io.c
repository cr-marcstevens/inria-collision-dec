#include <stdio.h>
#include <m4ri/m4ri.h>
#include "io.h"
#include "libisd.h"

int read_bit(FILE* fd) {
	int c;
	do {
		c = getc(fd);
	} while ((c != '0') && (c != '1') && (c != EOF));
	return c;
}

void CSD_from_file(mzd_t** HzeroT, unsigned int* w, unsigned int* N, word*** synds, FILE* fd) {
	unsigned int i, j;
	unsigned int n, r;
	int b;
	if(fscanf(fd, "%d %d %d %d", &n, &r, w, N) != 4) {
		fprintf(stderr, "Invalid instance file : could not read parameters\n");
		exit(2);
	}

	*synds = (word**) malloc(*N*sizeof(word*));
	for (i = 0; i < *N; ++i) {
		(*synds)[i] = (word*) calloc(bit_in_words(r), sizeof(word));
		for (j = 0; j < r; ++j) {
			b = read_bit(fd);
			if (b == EOF) {
				fprintf(stderr, "Invalid instance file : could not read syndromes");
				mzd_free(*HzeroT);
				fclose(fd);
				exit(2);
			}
			if (b == '1') {
				vec_set_bit((*synds)[i], j);
			}
		}
	}
	
	*HzeroT = mzd_init(n, r);
	for (i = 0; i < n; ++i) {
		for (j = 0; j < r; ++j) {
			b = read_bit(fd);
			if (b == EOF) {
				fprintf(stderr, "Invalid instance file : could not read matrix");
				mzd_free(*HzeroT);
				fclose(fd);
				exit(2);
			}
			mzd_write_bit(*HzeroT, i, j, b - '0');
		}
	}
}

void CSD_free(mzd_t* HzeroT, unsigned int N, word** synds) {
	unsigned int i;
	for (i = 0; i < N; ++i) {
		free(synds[i]);
	}
	free(synds);
	mzd_free(HzeroT);
}
