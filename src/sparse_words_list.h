#ifndef SPARSE_WORDS_LIST_H
#define SPARSE_WORDS_LIST_H

typedef struct sw {
	unsigned int weight;      // Weight of the word
	unsigned int *pos;        // Indexes of the non-null position of the word
	unsigned int synd_idx;    // Index of the target syndrom considered
	unsigned int synd_weight; // Weight of the syndrom of the word added to the target syndrom
	unsigned int sorted;      // 1 if pos is sorted, 0 instead

	struct sw *next;
} sw;

typedef sw sw_list;

sw_list* sw_new(unsigned int p);
sw* sw_filled_new(unsigned int synd_idx, unsigned int synd_weight, unsigned int p, ...);
void sw_sort(sw* h);
void sw_print(sw* h);
int sw_cmp(sw_list* h1, sw_list* h2);
void sw_free(sw_list* h);

void sw_list_add_array(sw_list** h, unsigned int synd_idx, unsigned int synd_weight, unsigned int p, unsigned short* columns);
void sw_list_append(sw_list** h, sw* new);
void sw_list_uniq(sw_list** h);
void sw_list_print(sw_list* h);
int sw_list_cmp(sw_list* h1, sw_list* h2);
void sw_list_free(sw_list* h);

#endif
