/*
 * BJMMtools.c
 *
 *  Created on: 27 juin 2013
 *      Author: Mathieu ARIA
 */

#include "m4ri/m4ri.h"
#include "BJMMtools.h"
#include <limits.h>

// see BJMMtools.h for method description

int next(word* L,unsigned short* indice,word* sums,unsigned int csize,unsigned int L_len, unsigned int w,int method){
	unsigned short pointer = w-1;
	unsigned int i,ii;
	switch (method){
	case 0:
			while(indice[pointer] == ((L_len/2)-w+pointer)){
				pointer--;
				if (pointer ==  USHRT_MAX) {
					return 0; // List fully build
				}
			}
			indice[pointer]+=1;
			for (i=(pointer+1); i<w; i++){
				indice[i]=(indice[i-1]+1);
			}
			// recomputing of partial sums
				if (pointer == 0) {
					for (i=0;i<csize;i++){
						sums[i*w]=L[indice[0]+i*L_len];
						for (ii=pointer+1;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
				else {
					for (i=0;i<csize;i++){
						for (ii=pointer;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
			return 1;
	case 1:
			while(indice[pointer] == (L_len+(-w+pointer)*2)){
				pointer--;
				if (pointer == USHRT_MAX) {
					return 0; // List fully build
				}
			}
			indice[pointer]+=2;
			for (i=(pointer+1); i<w; i++){
				indice[i]=(indice[i-1]+2);
			}
			// recomputing of partial sums
				if (pointer == 0) {
					for (i=0;i<csize;i++){
						sums[i*w]=L[indice[0]+i*L_len];
						for (ii=pointer+1;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
				else {
					for (i=0;i<csize;i++){
						for (ii=pointer;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
			return 1;
	case 2:
			while(indice[pointer] == ((L_len*3)/4-w+pointer)){
				pointer--;
				if (pointer == USHRT_MAX) {
					return 0; // List fully build
				}
			}
			if (indice[pointer]==((L_len/4)-1)){
				indice[pointer]=(L_len/2);
			}
			else{
				indice[pointer]+=1;
			}
			for (i=(pointer+1); i<w; i++){
				if (indice[i-1]==((L_len/4)-1)){
					indice[i]=(L_len/2);
				}
				else{
					indice[i]=indice[i-1]+1;
				}
			}
			// recomputing of partial sums
				if (pointer == 0) {
					for (i=0;i<csize;i++){
						sums[i*w]=L[indice[0]+i*L_len];
						for (ii=pointer+1;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
				else {
					for (i=0;i<csize;i++){
						for (ii=pointer;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
	case 3:
			while(indice[pointer] == (L_len-w+pointer)){
				pointer--;
				if (pointer == USHRT_MAX) {
					return 0; // List fully build
				}
			}
			if (indice[pointer]==((L_len/4)-1)){
				indice[pointer]=((L_len*3)/4);
			}
			else{
				indice[pointer]+=1;
			}
			for (i=(pointer+1); i<w; i++){
				if (indice[i-1]==((L_len/4)-1)){
					indice[i]=((L_len*3)/4);
				}
				else{
					indice[i]=indice[i-1]+1;
				}
			}
			// recomputing of partial sums
				if (pointer == 0) {
					for (i=0;i<csize;i++){
						sums[i*w]=L[indice[0]+i*L_len];
						for (ii=pointer+1;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
				else {
					for (i=0;i<csize;i++){
						for (ii=pointer;ii<w;ii++){
							sums[ii+i*w]=sums[(ii-1)+i*w]^L[indice[ii]+i*L_len];
						}
					}
				}
	default: printf("unknown building method");
			return 0;
	}
}

void h1store(S* table,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	word index = sums[w-1];
	if (table[index].indice==NULL){
		table[index].indice = malloc(w*sizeof(short));
		memcpy(table[index].indice,indice,w*sizeof(short));
		table[index].sum = malloc(csize*sizeof(word));
		for (i=0; i<csize; i++){
			memcpy(table[index].sum+i,sums+(w*(1+i)-1),sizeof(word));
		}
	}
	else{
		S* current = &table[index];
		while((*current).next!=NULL){
			current = (*current).next;
		}
		(*current).next= calloc(1,sizeof(S));
		(*((*current).next)).indice = malloc(w*sizeof(short));
		memcpy((*((*current).next)).indice,indice,w*sizeof(short));
		(*((*current).next)).sum = malloc(csize*sizeof(word));
		for (i=0; i<csize; i++){
			memcpy((*((*current).next)).sum+i,sums+(w*(1+i)-1),sizeof(word));
		}
	}
}

void freelist(S draw){
	if(draw.next != NULL){
		freelist(*(draw.next));
		free(draw.next);
	}
	free(draw.sum);
	free(draw.indice);
}

