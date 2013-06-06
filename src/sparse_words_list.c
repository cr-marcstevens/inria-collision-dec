#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "sparse_words_list.h"

sw_list* sw_list_new(sw_list* h, unsigned int p) {
	sw_list* new = (sw_list*) malloc(sizeof(sw_list));
	new->weight = p;
	new->pos = (unsigned int*) malloc(p*sizeof(unsigned int));
	new->synd_weight = 0;
	new->next = h;
	return new;
}

sw_list* sw_list_add(sw_list* h, unsigned int synd_idx, unsigned int synd_weight, unsigned int p, ...) {
	unsigned int i;
	va_list columns;

	sw_list* new = sw_list_new(h, p);

	va_start(columns, p);
	for (i = 0; i < p; ++i) {
		new->pos[i] = va_arg(columns, int);
	}
	va_end(columns);

	new->synd_idx = synd_idx;
	new->synd_weight = synd_weight;

	return new;
}

void sw_list_sort(sw_list* h) {
	int cmp(const void* a, const void* b) {
		return (*(int*) a > *(int*) b);
	}
	qsort(h->pos, h->weight, sizeof(int), cmp);
}

void sw_list_print(sw_list* h) {
	unsigned int i = 0;
	while(h) {
		printf("weight : %3d, synd_idx : %2d, synd_weight : %3d; ", h->weight, h->synd_idx, h->synd_weight);
		for (i = 0; i < h->weight; ++i) {
			printf("%5d ", h->pos[i]);
		}
		printf("\n");
		h = h->next;
	}
}

void sw_list_free(sw_list* h) {
	sw_list* next;
	while(h) {
		next = h->next;
		free(h->pos);
		free(h);
		h = next;
	}
}
