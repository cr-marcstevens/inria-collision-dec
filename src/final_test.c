#include <stdlib.h>
#include <stdarg.h>
#include "final_test.h"
#include "libisd.h"
#include <m4ri/m4ri.h>

static unsigned int r;
static unsigned int w;
static mzd_t* BT;
static mzd_t* Usecondmod;
static word** syndszero;
static word* eprimeBT;

static unsigned int eff_word_len;

void final_test_init(unsigned int local_r, unsigned int local_w, mzd_t* local_BT, mzd_t* local_Usecondmod, word** local_syndszero) {
	r = local_r;
	w = local_w;
	BT = local_BT;
	Usecondmod = local_Usecondmod;
	syndszero = local_syndszero;
	eff_word_len = min((unsigned int) BT->ncols, word_len);

	eprimeBT = (word*) malloc(bit_in_words(r) * sizeof(word));
}

/* Take a p-uplet of columns giving a zero symdrom on H' as an input and check if the weight of the coresponding syndrom on BT is <= w-p 
 * Returns -1 if weight > w-p
 * Else returns weight
 * */
int final_test(unsigned int synd_idx, unsigned int weight, unsigned int p, ... ) {
	va_list columns;
	int col;

	unsigned int i, j;

	/* generating e'*BT adding p lines of BT */
	va_start(columns, p);
	for (i = 0; i < bit_in_words(r); ++i) {
		eprimeBT[i] = syndszero[synd_idx][i];
	}
	for (j = 0; j < p; ++j) {
		col = va_arg (columns, int);
		for (i = 0; i < bit_in_words(r); ++i) {
			eprimeBT[i] ^= BT->rows[col][i];
		}
	}
	
	va_end(columns);

	/* applying gaussian elimination transformation to e'*BT (except eff_word_len first bits since we know their weight), and computing the weight of the result */
	for (i = 0; i < r - eff_word_len; ++i) {
		if (pscal(eprimeBT, Usecondmod->rows[i], r)) {
			++weight;
		}
		if (weight > w-p) {
			return -1;
		}
	}
	
	return weight;
}

int final_test_array(unsigned int synd_idx, unsigned int weight, unsigned int p, unsigned short* columns ) {

	unsigned int i, j;

	/* generating e'*BT adding p lines of BT */
	for (i = 0; i < bit_in_words(r); ++i) {
		eprimeBT[i] = syndszero[synd_idx][i];
	}
	for (j = 0; j < p; ++j) {
		for (i = 0; i < bit_in_words(r); ++i) {
			eprimeBT[i] ^= BT->rows[columns[j]][i];
		}
	}

	/* applying gaussian elimination transformation to e'*BT (except eff_word_len first bits since we know their weight), and computing the weight of the result */
	for (i = 0; i < r - eff_word_len; ++i) {
		if (pscal(eprimeBT, Usecondmod->rows[i], r)) {
			++weight;
		}
		if (weight > w-p) {
			return -1;
		}
	}

	return weight;
}

void final_test_free() {
	free(eprimeBT);
}
