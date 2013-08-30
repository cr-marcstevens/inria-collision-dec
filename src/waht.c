#include <stdlib.h>
#include "waht.h"
#include "syscall_macros.h"

waht waht_init(unsigned long long size) {
	return (waht) CALLOC(size, sizeof(waht_list));
}

void waht_store(waht L, word index, waht_elt value) {
	if (L[index] == NULL) {
		L[index] = (waht_list) malloc(2*sizeof(waht_elt));
		L[index][0] = 1;
		L[index][1] = value;
	}
	else {
		++L[index][0];
		L[index] = (waht_list) realloc(L[index], (1+L[index][0])*sizeof(waht_elt));
		L[index][L[index][0]] = value;
	}
}

waht_list waht_get(waht L, word index) {
	return L[index];
}

void waht_reset(waht L, unsigned long long size) {
	unsigned int i;
	for (i = 0; i < size; ++i) {
		free(L[i]);
		L[i] = NULL;
	}
}

void waht_free(waht L, int size) {
	waht_reset(L, size);
	free(L);
}
