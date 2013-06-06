#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "htable.h"

/*
 *	This structure stores integers indexed by integers
 *	Storing element a in index i is done this way :
 *	 - if htable[i] is empty, it receives a
 *	 - else if htable_col[htable[i]] is empty, it receives a
 *	 - else if htable_col[htable_col[htable[i]]] is empty, it receives a
 *	 - and so on
 */

static unsigned long long htable_size;
static int* htable = NULL;
static unsigned long long htable_col_size;
static int* htable_col = NULL;

void htable_init(unsigned long long max_idx, unsigned long long nb_elt) {
	htable_size = max_idx;
	htable_col_size = nb_elt;
	
	htable = (int*) malloc(htable_size * sizeof(int));
	htable_col = (int*) malloc(htable_col_size * sizeof(int));
}

void htable_store(unsigned int index, int value) {
	int current;
	if (htable[index] == NONE) {
		htable[index] = value;
		return;
	}
	current = htable[index];
	while (htable_col[current] != NONE){
		current = htable_col[current];
	}
	htable_col[current] = value;
}

int htable_get(unsigned int index) {
	return htable[index];
}

int htable_next(unsigned int index, int value) {
	(void) index;
	return htable_col[value];
}

void htable_reset() {
//	memset(d,0,10*sizeof(*d));
	memset(htable, NONE, htable_size*sizeof(*htable));
	memset(htable_col, NONE, htable_col_size*sizeof(*htable_col));
}

void htable_free() {
	free(htable);
	free(htable_col);
}

#define MAX_COUNTER 20
void htable_stats() {
	unsigned int i;
	int nb_htable = 0;
	int collision[MAX_COUNTER] = {0};
	int other = 0;
	int counter;
	int current;
	int max_counter = 0;
	int max_index = 0;

	for (i = 0; i < htable_size; ++i) {
		counter = 0;
		if (htable[i] != -1) {
			++nb_htable;
			counter = 1;
			current = htable[i];
			while(htable_col[current] != -1){
				current = htable_col[current];
				++counter;
			}
		}
		if (counter < MAX_COUNTER) {
			++collision[counter];
		}
		else {
			++other;
		}
		if (counter > max_counter) {
			max_counter = counter;
			max_index = i;
		}
		if (i == 0) {
			printf("Number of elts at idx 0 : %d\n", counter);
		}
	}
	printf("max index %d reached %d times\n", max_index, max_counter);

	int nb_htable_col = 0;
	for (i = 0; i < htable_col_size; ++i) {
		if (htable_col[i] != -1) {
			nb_htable_col++;
		}
	}

	printf("htable occupation ratio : %d/%lld : %3d%%\n", nb_htable, htable_size, (int) (100.0*nb_htable/htable_size));
	printf("htable_col occupation ratio : %d/%lld : %3d%%\n", nb_htable_col, htable_col_size, (int) (100.0*nb_htable_col/htable_col_size));
	printf("memory efficiency : %3d/%6lld : %3lld%%\n", nb_htable+nb_htable_col, htable_size+htable_col_size, (nb_htable+nb_htable_col)/(htable_size+htable_col_size));
	printf("collisions : \n");
	for (i = 0; i < MAX_COUNTER; ++i) {
		printf("%d : %d\n", i, collision[i]);
	}
	printf("+ : %d\n", other);
}
