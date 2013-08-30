#include <stdlib.h>
#include "counterht_nocol.h"
#include "syscall_macros.h"

counterht counterht_init(unsigned long long size) {
	return (counterht) CALLOC(size, sizeof(counter));
}

void counterht_store(counterht L, word index, counter value) {
	L[index] = value;
}

counter counterht_get(counterht L, word index) {
	return L[index];
}

counter counterht_next(counterht L, word index, counter c) {
	(void) L;
	(void) index;
	(void) c;
	return NONE;
}

void counterht_reset(counterht L, unsigned long long L_size) {
	unsigned int i;
	for (i = 0; i < L_size; ++i) {
		L[i] = NONE;
	}
}

void counterht_free(counterht L) {
	free(L);
}
