#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "sparse_words_list.h"

sw* sw_new(unsigned int p) {
	sw_list* new = (sw_list*) malloc(sizeof(sw_list));
	new->weight = p;
	new->pos = (unsigned int*) malloc(p*sizeof(unsigned int));
	new->synd_weight = 0;
	new->sorted = 0;
	new->next = NULL;
	return new;
}

sw* sw_filled_new(unsigned int synd_idx, unsigned int synd_weight, unsigned int p, ...) {
	unsigned int i;
	va_list columns;

	sw* new = sw_new(p);

	va_start(columns, p);
	for (i = 0; i < p; ++i) {
		new->pos[i] = va_arg(columns, int);
	}
	va_end(columns);

	new->synd_idx = synd_idx;
	new->synd_weight = synd_weight;

	return new;
}

void sw_list_append(sw_list** h, sw* new) {
	new->next = *h;
	*h = new;
}

void sw_list_add_array(sw_list** h, unsigned int synd_idx, unsigned int synd_weight, unsigned int p, unsigned short* columns) {
	unsigned int i;
	sw* new = sw_new(p);
	for (i = 0; i < p; ++i) {
		new->pos[i] = *(columns+i);
	}
	new->synd_idx = synd_idx;
	new->synd_weight = synd_weight;
	sw_list_append(h, new);
}

void sw_list_sort(sw_list* h) {
	int cmp(const void* a, const void* b) {
		return (*(int*) a > *(int*) b);
	}
	qsort(h->pos, h->weight, sizeof(int), cmp);
}

void sw_print(sw* w) {
	unsigned int i = 0;
	printf("weight : %3d, synd_idx : %2d, synd_weight : %3d; ", w->weight, w->synd_idx, w->synd_weight);
	for (i = 0; i < w->weight; ++i) {
		printf("%5d ", w->pos[i]);
	}
	printf("\n");
}

void sw_list_print(sw_list* h) {
	while(h) {
		sw_print(h);
		h = h->next;
	}
}

void sw_sort(sw* w) {
	if (w->sorted) {
		return;
	}
	int cmp(const void* a, const void* b) {
		return *(unsigned int*)a-*(unsigned int*)b;
	}
	qsort(w->pos, w->weight, sizeof(unsigned int), cmp);
	w->sorted = 1;
}

int sw_cmp(sw* w1, sw* w2) {
	if (w1->weight != w2->weight) {
		return -1;
	}
	else {
		sw_sort(w1);
		sw_sort(w2);
		unsigned int i;
		for (i = 0; i < w1->weight; ++i) {
			if (w1->pos[i] != w2->pos[i]) {
				return -1;
			}
		}
		return 0;
	}
}

int belongs(sw* w, sw_list* h) {
	while (h) {
		if (sw_cmp(w, h) == 0) {
			return 1;
		}
		h = h->next;
	}
	return 0;
}

/* Filter the list removing multiples */
void sw_list_uniq(sw_list** h) {
	sw_list* ptr = *h;
	sw_list* g = NULL;
	sw_list* next;
	while (ptr) {
		next = ptr->next;
		if (!belongs(ptr, g)) {
			sw_list_append(&g, ptr);
		}
		else {
			sw_free(ptr);
		}
		ptr = next;
	}
	*h = g;
}

void sw_free(sw* w) {
	free(w->pos);
	free(w);
}

void sw_list_free(sw_list* h) {
	sw_list* next;
	while(h) {
		next = h->next;
		sw_free(h);
		h = next;
	}
}
