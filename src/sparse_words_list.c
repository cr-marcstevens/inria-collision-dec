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
		new->pos[i] = va_arg(columns, unsigned int);
	}
	va_end(columns);

	new->synd_idx = synd_idx;
	new->synd_weight = synd_weight;

	return new;
}

sw* sw_filled_new_array(unsigned int synd_idx, unsigned int synd_weight, unsigned int p, unsigned short* columns) {
	unsigned int i;
	sw* new = sw_new(p);

	for (i = 0; i < p; ++i) {
		new->pos[i] = columns[i];
	}

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

unsigned int sw_list_len(sw_list* h) {
	unsigned int count = 0;
	while (h) {
		++count;
		h = h->next;
	}
	return count;
}

/*
sw_list* sw_list_sort(sw_list* h) {
	printf("SORTING\n");
	if (h == NULL || h->next == NULL) {
		return h;
	}
	sw_list* head = NULL;
	sw_list* cur;
	sw_list* p;
	while (h) {
		printf("HEAD : \n");
		sw_list_print(head);
		printf("\n");
		
		cur = h;
		printf("Checking\n");
		sw_print(cur);
		h = h->next;
		if (head == NULL || sw_cmp(cur, head) < 0) {
			if (head == NULL) {
				printf("head is NULL\n");
			}
			else {
				sw_print(cur);
				sw_print(head);
			}
			printf("Inserting at head\n");
			cur->next = head;
			head = cur;
		}
		else {
			printf("Inserting in the list\n");
			p = head;
			while (p) {
				printf("comparing\n");
				printf("cur : "); sw_print(cur);
				printf("p   : "); sw_print(p);
				printf("res : %d\n", sw_cmp(cur, p));
				if (p->next == NULL || sw_cmp(cur, p) < 0) {
					if (p->next == NULL) {
						printf("Inserting at end\n");
					}
					else {
						printf("Inserting middle\n");
					}
					cur->next = p->next;
					p->next = cur;
					break;
				}
				p = p->next;
			}
		}
	}
	printf("END SORTING\n");
	return head;
}
*/

/* http://www.geeksforgeeks.org/merge-sort-for-linked-list/ */
void sw_list_split(sw_list* h, sw_list** h1, sw_list** h2) {
  sw_list* fast;
  sw_list* slow;
  if (h==NULL || h->next==NULL)
  {
    /* length < 2 */
    *h1 = h;
    *h2 = NULL;
  }
  else
  {
    slow = h;
    fast = h->next;
 
    /* Advance 'fast' two nodes, and advance 'slow' one node */
    while (fast != NULL)
    {
      fast = fast->next;
      if (fast != NULL)
      {
        slow = slow->next;
        fast = fast->next;
      }
    }
 
    /* 'slow' is before the midpoint in the list, so split it in two
      at that point. */
    *h1 = h;
    *h2 = slow->next;
    slow->next = NULL;
  }
}

sw_list* sw_list_merge(sw_list* h1, sw_list* h2) {
	  sw_list* res = NULL;
 
  /* Base cases */
  if (h1 == NULL)
     return(h2);
  else if (h2==NULL)
     return(h1);
 
  /* Pick either a or b, and recur */
  if (sw_cmp(h1, h2) < 0)
  {
     res = h1;
     res->next = sw_list_merge(h1->next, h2);
  }
  else
  {
     res = h2;
     res->next = sw_list_merge(h1, h2->next);
  }
  return(res);
}

/* This may be not the most efficient because we have to walk through the list to split it */
void sw_list_sort(sw_list** h) {
  sw_list* head = *h;
  sw_list* a;
  sw_list* b;
 
  /* Base case -- length 0 or 1 */
  if ((head == NULL) || (head->next == NULL))
  {
    return;
  }
 
  /* Split head into 'a' and 'b' sublists */
  sw_list_split(head, &a, &b); 
 
  /* Recursively sort the sublists */
  sw_list_sort(&a);
  sw_list_sort(&b);
 
  /* answer = merge the two sorted lists together */
  *h = sw_list_merge(a, b);
}

/*
fonction trier(p, n)
    Q := n/2 (division entière)
    P := n-Q
    si P >= 2
        q := trier(p, P)
        si Q >= 2 trier(q, Q)
    sinon
        q := p.suivant
    fin
    q := fusionner(p, P, q, Q)
    renvoyer q
fin
 
fonction fusionner(p, P, q, Q)
    répéter indéfiniment
        si valeur(p.suivant) > valeur(q.suivant) 
            déplacer le maillon q.suivant après le maillon p
            si Q = 1 quitter la boucle
            Q := Q-1
        sinon
            si P = 1
                tant que Q >= 1
                    q := q.suivant
                    Q := Q-1
                fin
                quitter la boucle
            fin
            P := P-1
        fin
        p := p.suivant
    fin
    renvoyer q
fin
*/

int sw_list_sorted(sw_list* h) {
	if (h == NULL) {
		return 1;
	}
	while (h->next) {
		if (sw_cmp(h, h->next) > 0) {
			return 0;
		}
		h = h->next;
	}
	return 1;
}

void sw_print(sw* w) {
	if (w == NULL) {
		printf("NULL\n");
		return;
	}
	unsigned int i = 0;
	printf("weight : %3d, synd_idx : %2d, synd_weight : %3d; ", w->weight, w->synd_idx, w->synd_weight);
	for (i = 0; i < w->weight; ++i) {
		printf("%5d ", w->pos[i]);
	}
	printf("\n");
}

void sw_list_print(sw_list* h) {
	if (h == NULL) {
		printf("NULL\n");
		return;
	}
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
			if (w1->pos[i] < w2->pos[i]) {
				return -1;
			}
			else if(w1->pos[i] > w2->pos[i]) {
				return 1;
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
int sw_list_uniq(sw_list** h) {
	sw_list* ptr = *h;
	sw_list* g = NULL;
	sw_list* next;
	int counter = 0;
	while (ptr) {
		next = ptr->next;
		if (!belongs(ptr, g)) {
			sw_list_append(&g, ptr);
			++counter;
		}
		else {
			sw_free(ptr);
		}
		ptr = next;
	}
	*h = g;
	return counter;
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
