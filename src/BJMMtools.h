/*
 * BJMMtools.h
 *
 *  Created on: 27 juin 2013
 *      Author: Mathieu Aria
 */

typedef struct draw {
	int* indice;
	struct draw* next;
} S ;

int next(short* indice,unsigned int L_len, unsigned int w,int method);
