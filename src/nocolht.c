#include <stdlib.h>
#include "nocolht.h"
#include "syscall_macros.h"

nocolht nocolht_init(unsigned long long size) {
	return (nocolht) CALLOC(size, sizeof(nocolht_elt));
}

void nocolht_store(nocolht L, word index, nocolht_elt value) {
	L[index] = value;
}

nocolht_elt nocolht_get(nocolht L, word index) {
	return L[index];
}

void nocolht_reset(nocolht L, unsigned long long L_size) {
	unsigned int i;
	for (i = 0; i < L_size; ++i) {
		L[i] = 0; /* TODO : 0 is a value that could be stored in the table */
	}
}

void nocolht_free(nocolht L) {
	free(L);
}
