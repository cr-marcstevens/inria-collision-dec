#ifndef COUNTERHT_COL2_H
#define COUNTERHT_COL2_H
#include "counterht_col2.h"
#include "syscall_macros.h"

counter counter_container_open(counter_container* c) {
	return *c;
}

counterht counterht_init(unsigned long long nb_of_hash_val, unsigned long long nb_of_insert) {
	counterht L = (counterht) malloc(sizeof(*L));
	L->main_table = (counter*) MALLOC(nb_of_hash_val*sizeof(counter));
	L->col_table = (counter*) MALLOC(nb_of_insert*sizeof(counter));
	L->main_table_size = nb_of_hash_val;
	L->col_table_size = nb_of_insert;
	return L;
}

void counterht_store(counterht L, word index, counter value) {
	counter_container* maint = L->main_table;
	if (maint[index] == NONE) {
		maint[index] = value;
		return;
	}
	counter current = maint[index];
	counter* colt = L->col_table;
	while (colt[current] != NONE){
		current = colt[current];
	}
	colt[current] = value;
}

counter_container* counterht_get(counterht L, word index) {
	if (L->main_table[index] == NONE) {
		return NULL;
	}
	return &L->main_table[index];
}

counter_container* counterht_next(counterht L, word index, counter_container* c) {
	(void) index;
	if (L->col_table[counter_container_open(c)] == NONE) {
		return NULL;
	}
	return &L->col_table[counter_container_open(c)];
}

void counterht_reset(counterht L, unsigned long long L_size) {
	(void) L_size;
	memset(L->main_table, NONE, L->main_table_size*sizeof(counter));
	memset(L->col_table, NONE, L->col_table_size*sizeof(counter));
}

void counterht_free(counterht L) {
	free(L->main_table);
	free(L->col_table);
	free(L);
}
#endif
