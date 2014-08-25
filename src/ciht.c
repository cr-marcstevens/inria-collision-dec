#include <stdlib.h>
#include "ciht.h"
#include "syscall_macros.h"

ciht_t ciht_init(unsigned long long nb_of_hash_val, unsigned int p, unsigned long long nb_of_insert) {
	(void) nb_of_insert;
	ciht_t ret = (ciht_t) MALLOC(p*nb_of_hash_val*sizeof(ci_t));
	return ret;
}

void ciht_store(ciht_t L, word index, ci_t* value, unsigned int p) {
	memcpy(L+(p*index), value, p*sizeof(ci_t));
}

ci_t* ciht_get(ciht_t L, word index, unsigned int p) {
	if (L[p*index] == NONE) {
		return NULL;
	}
	return L+p*index;
}

void ciht_reset(ciht_t L, unsigned long long L_size, unsigned int p) {
	memset(L, NONE, p*L_size*sizeof(ci_t));
}

void ciht_free(ciht_t L) {
	free(L);
}
