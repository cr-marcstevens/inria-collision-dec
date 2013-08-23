#include <stdlib.h>
#include <stdarg.h>
#include "m4ri/m4ri.h"
#include "iaht.h"
#include "syscall_macros.h"

/*
 * iaht (integer array has table) is an array hash table containing list of integer x-uplet
 * First element of each lists is the number
 * of x-uplet in the list. Next x columns are the columns of the first
 * x-uplet and so on...
 */

iaht iaht_init(unsigned long long size) {
	return (iaht) CALLOC(size, sizeof(iaht_list));
}


void iaht_store(iaht L, word index, int x, ...) {
	va_list indices;
	int indice;
	int i;
	iaht_list current = L[index];
	va_start(indices, x);
	if (current == NULL) {
		current = (iaht_list) malloc((1+x)*sizeof(iaht_elt));
		current[0] = 1;
	}
	else {
		++current[0];
		current = (iaht_list) realloc(current, (x*current[0] + 1)*sizeof(iaht_elt));
	}
	for (i = 0; i < x; ++i) {
		indice = va_arg(indices, iaht_elt);
		current[x*(current[0]-1)+1+i] = indice;
	}
	L[index] = current;
}

iaht_list iaht_get(iaht L, word index) {
	return L[index];
}

void iaht_reset(iaht L, unsigned long long size) {
	unsigned int i;
	for (i = 0; i < size; ++i) {
		free(L[i]);
		L[i] = NULL;
	}
}

void iaht_free(iaht L, int size) {
	iaht_reset(L, size);
	free(L);
}
