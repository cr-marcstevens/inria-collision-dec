/*
 * BJMMtools.h
 *
 *  Created on: 27 juin 2013
 *      Author: Mathieu Aria
 */

typedef struct draw {
	word* sum;
	unsigned short* indice;
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

int next2(word* L,unsigned short* indice,word* sums,unsigned int csize,unsigned int L_len, unsigned int w,int method);

void h1store(S* table,int shift,unsigned short* indice,word* sums,unsigned int w,unsigned int csize);

void fusion1(S* EStep1,word target,S* table,int shift1,int shift2,int shift3,unsigned short* indice,word* sums,unsigned int w,unsigned int csize);

void tri(unsigned short* dest,unsigned short* s1,unsigned short* s2,unsigned short size1,unsigned short size2);

void freelist(S draw);
