#include <stdlib.h>
#include "counterht.h"
#include "syscall_macros.h"

counter counter_container_open(counter_container* c) {
	return *c;
}

counterht counterht_init(unsigned long long nb_of_hash_val, unsigned long long nb_of_insert) {
	(void) nb_of_insert;
	counterht ret = (counterht) CALLOC(nb_of_hash_val, sizeof(*ret));
	return ret;
}

void counterht_store(counterht L, word index, counter value) {
	L[index] = value;
}

counter_container* counterht_get(counterht L, word index) {
	if (L[index] == NONE) {
		return NULL;
	}
	return &L[index];
}

counter_container* counterht_next(counterht L, word index, counter_container* c) {
	(void) L;
	(void) index;
	(void) c;
	return NULL;
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
