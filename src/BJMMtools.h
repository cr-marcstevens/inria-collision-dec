/*
 * BJMMtools.h
 *
 *  Created on: 27 juin 2013
 *      Author: Mathieu Aria
 */

typedef struct draw {
	word* sum;
	short* indice;
	struct draw* next;
} S ;

/*
 *  About next:
 *  next can be called with 4 method (from 0 to 3) which share the matrix L with different patterns:
 *
 *  method 0 : [x------x||0------0] first half
 *	method 1 : [x0x0x0x0x0x0x0x0x0] odd number
 *	method 2 : [x-x||0-0||x-x||0-0]	first and third quarter
 *	method 3 : [x-x||0------0||x-x]	first and final quarter
 */

int next(word* L,unsigned short* indice,word* sums,unsigned int csize,unsigned int L_len, unsigned int w,int method);

void h1store(S* table,unsigned short* indice,word* sums,unsigned int w,unsigned int csize);

void freelist(S draw);
