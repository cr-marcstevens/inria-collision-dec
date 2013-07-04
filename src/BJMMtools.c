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
						for (ii=1;ii<w;ii++){
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
						for (ii=1;ii<w;ii++){
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
						for (ii=1;ii<w;ii++){
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
						for (ii=1;ii<w;ii++){
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

void h1store(S* table,int shift,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	word index = (sums[w-1] >> shift);
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

int next2(word* L,unsigned short* indice,word* sums,unsigned int csize,unsigned int L_len, unsigned int w,int method){
	unsigned short pointer = w-1;
	unsigned int i,ii;
	switch (method){
	case 0:
			while(indice[pointer] == ((L_len)-w+pointer)){
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
						for (ii=1;ii<w;ii++){
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
			while(indice[pointer] == (L_len+(-w+pointer)*2+1)){
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
						for (ii=1;ii<w;ii++){
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
			while(indice[pointer] == (L_len-w+pointer)){
				pointer--;
				if (pointer == USHRT_MAX) {
					return 0; // List fully build
				}
			}
			if (indice[pointer]==((L_len/2)-1)){
				indice[pointer]=((L_len*3)/4);
			}
			else{
				indice[pointer]+=1;
			}
			for (i=(pointer+1); i<w; i++){
				if (indice[i-1]==((L_len/2)-1)){
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
						for (ii=1;ii<w;ii++){
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
			while(indice[pointer] == (((L_len*3)/4)-w+pointer)){
				pointer--;
				if (pointer == USHRT_MAX) {
						return 0; // List fully build
				}
			}
			indice[pointer]+=1;
			for (i=(pointer+1); i<w; i++){
					indice[i]=indice[i-1]+1;
			}
			// recomputing of partial sums
				if (pointer == 0) {
					for (i=0;i<csize;i++){
						sums[i*w]=L[indice[0]+i*L_len];
						for (ii=1;ii<w;ii++){
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

void fusion1(S* EStep1,word target,S* table,int shift1,int shift2,int shift3,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	word sumr1;
	word index = ((sums[w-1]>>shift1)^target);
	if (table[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		sumr1 = (((sums[w-1]^(table[index].sum[0]))<<shift2)>>shift3);
		if(EStep1[sumr1].indice==NULL){			//there is not yet a solution at sumr1
			EStep1[sumr1].sum = malloc(csize*sizeof(word));
			EStep1[sumr1].indice = malloc(w*sizeof(short));
			for (i=0; i<csize; i++){
				*(EStep1[sumr1].sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
			}
			tri((*EStep1).indice,table[index].indice,indice,w/2,w/2);
		}
		else{			//there are already one or more solutions at sumr1
			S* current1 = &EStep1[sumr1];
			while((*current1).next!=NULL){
				current1 = (*current1).next;
			}
			(*current1).next= calloc(1,sizeof(S));
			(*((*current1).next)).indice = malloc(w*sizeof(short));
			tri((*((*current1).next)).indice,table[index].indice,indice,w/2,w/2);
			(*((*current1).next)).sum = malloc(csize*sizeof(word));
			for (i=0; i<csize; i++){
				*((*((*current1).next)).sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
			}
		}
		S* current2 = &table[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			sumr1 = (((sums[w-1]^((*current2).sum[0]))<<shift2)>>shift3);

			if(EStep1[sumr1].indice==NULL){			//there is not yet a solution at sumr1
				EStep1[sumr1].sum = malloc(csize*sizeof(word));
				EStep1[sumr1].indice = malloc(w*sizeof(short));
				for (i=0; i<csize; i++){
					*(EStep1[sumr1].sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
				}
				tri((*EStep1).indice,(*current2).indice,indice,w/2,w/2);
			}
			else{			//there are already one or more solutions at sumr1
				S* current1 = &EStep1[sumr1];
				while((*current1).next!=NULL){
					current1 = (*current1).next;
				}
				(*current1).next= calloc(1,sizeof(S));
				(*((*current1).next)).indice = malloc(w*sizeof(short));
				tri((*((*current1).next)).indice,(*current2).indice,indice,w/2,w/2);
				(*((*current1).next)).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*((*current1).next)).sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
				}
			}
		}
	}
}

void tri(unsigned short* dest,unsigned short* s1,unsigned short* s2,unsigned short size1,unsigned short size2){
	unsigned short p1 =0;
	unsigned short p2 =0;
	while (p1 < size1 && p2 < size2){
		if (*(s1+p1) < *(s2+p2)){
			*(dest+p1+p2)= *(s1+p1);
			p1++;
		}
		else{
			*(dest+p1+p2)= *(s2+p2);
			p2++;
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

