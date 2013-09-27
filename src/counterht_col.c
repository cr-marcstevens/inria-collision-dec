#include <stdlib.h>
#include "counterht_col.h"
#include "syscall_macros.h"

counter counter_container_open(counter_container* c) {
	return c->value;
}

counterht counterht_init(unsigned long long nb_of_hash_val, unsigned long long nb_of_insert) {
	(void) nb_of_insert;
	counterht ret = (counterht) CALLOC(nb_of_hash_val, sizeof(*ret));
	return ret;
}

void counterht_store(counterht L, word index, counter value) {
	counter_container* ccont = (counter_container*) malloc(sizeof(counter_container));
	ccont->value = value;
	ccont->next = L[index];
	L[index] = ccont;
}

counter_container* counterht_get(counterht L, word index) {
	return L[index];
}

counter_container* counterht_next(counterht L, word index, counter_container* c) {
	(void) L;
	(void) index;
	return c->next;
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
