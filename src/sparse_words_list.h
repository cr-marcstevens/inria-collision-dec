#ifndef SPARSE_WORDS_LIST_H
#define SPARSE_WORDS_LIST_H

typedef struct sw {
	unsigned int weight;      // Weight of the word
	unsigned int *pos;        // Indexes of the non-null position of the word
	unsigned int synd_idx;    // Index of the target syndrom considered
	unsigned int synd_weight; // Weight of the syndrom of the word added to the target syndrom

	struct sw *next;
} sw_list;

sw_list* sw_list_new(sw_list* h, unsigned int p);
sw_list* sw_list_add(sw_list* h, unsigned int synd_idx, unsigned int synd_weight, unsigned int p, ...);
sw_list* sw_list_add_array(sw_list* h, unsigned int synd_idx, unsigned int synd_weight, unsigned int p,unsigned short* columns);
void sw_list_sort(sw_list* h);
void sw_list_print(sw_list* h);
void sw_list_free(sw_list* h);

#endif
