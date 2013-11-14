/*
 * BJMMtools.c
 *
 *  Created on: 27 juin 2013
 *      Author: Mathieu ARIA
 */

#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <m4ri/m4ri.h>
#include "sparse_words_list.h"
#include "final_test.h"
#include "BJMMtools.h"
#include "libisd.h"
#include "measure.h"
#include "cpucycles/cpucycles.h"


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
			while(indice[pointer] == (L_len+(pointer-w)*2)){
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
			return 1;
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
			return 1;
	default: printf("unknown building method \n");
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
					//printf("\n"); // barre de suivi.
					return 0; // List fully build

				}
			}
			indice[pointer]+=1;
			for (i=(pointer+1); i<w; i++){
				indice[i]=(indice[i-1]+1);
			}
			// recomputing of partial sums
				if (pointer == 0) {
					//printf("."); // barre de suivi.
					fflush(stdout);
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
					//printf("\n"); // barre de suivi.
					return 0; // List fully build
				}
			}
			indice[pointer]+=2;
			for (i=(pointer+1); i<w; i++){
				indice[i]=(indice[i-1]+2);
			}
			// recomputing of partial sums
				if (pointer == 0) {
					//printf("."); // barre de suivi.
					fflush(stdout);
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
					//printf("\n"); // barre de suivi.
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
					//printf("."); // barre de suivi.
					fflush(stdout);
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

	case 3:
			while(indice[pointer] == (((L_len*3)/4)-w+pointer)){
				pointer--;
				if (pointer == USHRT_MAX) {
					//printf("\n"); // barre de suivi.
					return 0; // List fully build
				}
			}
			indice[pointer]+=1;
			for (i=(pointer+1); i<w; i++){
					indice[i]=indice[i-1]+1;
			}
			// recomputing of partial sums
				if (pointer == 0) {
					//printf("."); // barre de suivi.
					fflush(stdout);
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
	default: printf("warning: unknown building method \n");
			return 0;
			}
}

void fusionstore1(S* EStep1,word target,S* table,int shift1,int shift2,int shift3,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	word sumr1;
	word index = ((sums[w/2-1]>>shift1)^target);
	if (table[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		sumr1 = (((sums[w/2-1]^(table[index].sum[0]))<<shift2)>>shift3);
		if(EStep1[sumr1].indice==NULL){			//there is not yet a solution at sumr1
			EStep1[sumr1].sum = malloc(csize*sizeof(word));
			EStep1[sumr1].indice = malloc(w*sizeof(short));
			for (i=0; i<csize; i++){
				*(EStep1[sumr1].sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
			}
			Sort(EStep1[sumr1].indice,table[index].indice,indice,w/2,w/2);
		}
		else{			//there are already one or more solutions at sumr1
			S* current1 = &EStep1[sumr1];
			while((*current1).next!=NULL){
				current1 = (*current1).next;
			}
			(*current1).next= calloc(1,sizeof(S));
			(*((*current1).next)).indice = malloc(w*sizeof(short));
			Sort((*((*current1).next)).indice,table[index].indice,indice,w/2,w/2);
			(*((*current1).next)).sum = malloc(csize*sizeof(word));
			for (i=0; i<csize; i++){
				*((*((*current1).next)).sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
			}
		}

		S* current2 = &table[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			sumr1 = (((sums[w/2-1]^((*current2).sum[0]))<<shift2)>>shift3);

			if(EStep1[sumr1].indice==NULL){			//there is not yet a solution at sumr1
				EStep1[sumr1].sum = malloc(csize*sizeof(word));
				EStep1[sumr1].indice = malloc(w*sizeof(short));
				for (i=0; i<csize; i++){
					*(EStep1[sumr1].sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
				}
				Sort(EStep1[sumr1].indice,(*current2).indice,indice,w/2,w/2);
			}
			else{			//there are already one or more solutions at sumr1
				S* current1 = &EStep1[sumr1];
				while((*current1).next!=NULL){
					current1 = (*current1).next;
				}
				(*current1).next= calloc(1,sizeof(S));
				(*((*current1).next)).indice = malloc(w*sizeof(short));
				Sort((*((*current1).next)).indice,(*current2).indice,indice,w/2,w/2);
				(*((*current1).next)).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*((*current1).next)).sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
				}
			}

		}
	}
}

void fusiongive1(S* answer,word target,S* table,int shift1,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	word index = ((sums[w/2-1]>>shift1)^target);
	if (table[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		(*answer).sum = malloc(csize*sizeof(word));
		(*answer).indice = malloc(w*sizeof(short));
		for (i=0; i<csize; i++){
			*((*answer).sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
		}
		Sort((*answer).indice,table[index].indice,indice,w/2,w/2);
		S* current1 = answer;
		S* current2 = &table[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			(*current1).next=calloc(1,sizeof(S));
			current1 = (*current1).next;
			(*current1).sum = malloc(csize*sizeof(word));
			(*current1).indice = malloc(w*sizeof(short));
			for (i=0; i<csize; i++){
				*((*current1).sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
			}
			Sort((*current1).indice,(*current2).indice,indice,w/2,w/2);
		}
	}
}

void FusionFilterStore64(S* AnswerList, S* StockedE,S* OnTheFlyE,word target,int shift1,int shift2,int eff_word_len,unsigned int l,unsigned int l2,unsigned int l3,unsigned int w,unsigned int w2,unsigned int csize){
	unsigned int i;
	word sumr;
	word index = (((((*OnTheFlyE).sum[0])<<shift1)>>shift2)^target); // OnTheFlyE first draw
	if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		sumr = ((((*OnTheFlyE).sum[0]^(StockedE[index].sum[0]))<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3)));
		if(AnswerList[sumr].indice==NULL){			//there is not yet a solution at sumr1
			AnswerList[sumr].indice = malloc(w2*sizeof(short));
			if(SortFilter(AnswerList[sumr].indice,(*OnTheFlyE).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
				AnswerList[sumr].sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*(AnswerList[sumr].sum+i)=((*OnTheFlyE).sum[i]^(StockedE[index].sum[i]));
				}
			}
			else{	//filtered solution
				free(AnswerList[sumr].indice);
				AnswerList[sumr].indice = NULL;
			}
		}
		else{			//there are already one or more solutions at sumr1
			S* current1 = &AnswerList[sumr];
			while((*current1).next!=NULL){
				current1 = (*current1).next;
			}
			(*current1).next= calloc(1,sizeof(S));
			(*((*current1).next)).indice = malloc(w2*sizeof(short));
			if(SortFilter((*((*current1).next)).indice,(*OnTheFlyE).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
				(*((*current1).next)).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*((*current1).next)).sum+i)=((*OnTheFlyE).sum[i]^(StockedE[index].sum[i]));
				}
			}
			else{	//filtered solution
				free((*((*current1).next)).indice);
				free((*current1).next);
			}
		}

		S* current2 = &(StockedE[index]);
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			sumr = ((((*OnTheFlyE).sum[0]^((*current2).sum[0]))<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3)));

			if(AnswerList[sumr].indice==NULL){			//there is not yet a solution at sumr1
				AnswerList[sumr].indice = malloc(w2*sizeof(short));
				if(SortFilter(AnswerList[sumr].indice,(*OnTheFlyE).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
					AnswerList[sumr].sum = malloc(csize*sizeof(word));
					for (i=0; i<csize; i++){
						*(AnswerList[sumr].sum+i)=((*OnTheFlyE).sum[i]^((*current2).sum[i]));
					}
				}
				else{	//filtered solution
					free(AnswerList[sumr].indice);
					AnswerList[sumr].indice = NULL;
				}
			}
			else{			//there are already one or more solutions at sumr1
				S* current1 = &AnswerList[sumr];
				while((*current1).next!=NULL){
					current1 = (*current1).next;
				}
				(*current1).next= calloc(1,sizeof(S));
				(*((*current1).next)).indice = malloc(w2*sizeof(short));
				if(SortFilter((*((*current1).next)).indice,(*OnTheFlyE).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
					(*((*current1).next)).sum = malloc(csize*sizeof(word));
					for (i=0; i<csize; i++){
						*((*((*current1).next)).sum+i)=((*OnTheFlyE).sum[i]^((*current2).sum[i]));
					}
				}
				else{	//filtered solution
					free((*((*current1).next)).indice);
					free((*current1).next);
				}
			}
		}
	}
	S* current3 = OnTheFlyE;
	while((*current3).next != NULL){ // OnTheFlyE next draws
		current3 = (*current3).next;
		word index = (((((*current3).sum[0])<<shift1)>>shift2)^target);
		if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
			sumr = ((((*current3).sum[0]^(StockedE[index].sum[0]))<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3)));

			if(AnswerList[sumr].indice==NULL){			//there is not yet a solution at sumr1
				AnswerList[sumr].indice = malloc(w2*sizeof(short));
				if(SortFilter(AnswerList[sumr].indice,(*current3).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
					AnswerList[sumr].sum = malloc(csize*sizeof(word));
					for (i=0; i<csize; i++){
						*(AnswerList[sumr].sum+i)=((*current3).sum[i]^(StockedE[index].sum[i]));
					}
				}
				else{	//filtered solution
					free(AnswerList[sumr].indice);
					AnswerList[sumr].indice = NULL;
				}
			}
			else{			//there are already one or more solutions at sumr1
				S* current1 = &AnswerList[sumr];
				while((*current1).next!=NULL){
					current1 = (*current1).next;
				}
				(*current1).next= calloc(1,sizeof(S));
				(*((*current1).next)).indice = malloc(w2*sizeof(short));
				if(SortFilter((*((*current1).next)).indice,(*current3).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
					(*((*current1).next)).sum = malloc(csize*sizeof(word));
					for (i=0; i<csize; i++){
						*((*((*current1).next)).sum+i)=((*current3).sum[i]^(StockedE[index].sum[i]));
					}
				}
				else{	//filtered solution
					free((*((*current1).next)).indice);
					free((*current1).next);
				}
			}


			S* current2 = &StockedE[index];
			while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
				current2 = (*current2).next;
				sumr = ((((*current3).sum[0]^((*current2).sum[0]))<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3)));

				if(AnswerList[sumr].indice==NULL){			//there is not yet a solution at sumr1
					AnswerList[sumr].indice = malloc(w2*sizeof(short));
					if(SortFilter(AnswerList[sumr].indice,(*current3).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
						AnswerList[sumr].sum = malloc(csize*sizeof(word));
						for (i=0; i<csize; i++){
							*(AnswerList[sumr].sum+i)=((*current3).sum[i]^((*current2).sum[i]));
						}
					}
					else{	//filtered solution
						free(AnswerList[sumr].indice);
						AnswerList[sumr].indice = NULL;
					}
				}
				else{			//there are already one or more solutions at sumr1
					S* current1 = &AnswerList[sumr];
					while((*current1).next!=NULL){
						current1 = (*current1).next;
					}
					(*current1).next= calloc(1,sizeof(S));
					(*((*current1).next)).indice = malloc(w2*sizeof(short));
					if(SortFilter((*((*current1).next)).indice,(*current3).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
						(*((*current1).next)).sum = malloc(csize*sizeof(word));
						for (i=0; i<csize; i++){
							*((*((*current1).next)).sum+i)=((*current3).sum[i]^((*current2).sum[i]));
						}
					}
					else{	//filtered solution
						free((*((*current1).next)).indice);
						free((*current1).next);
					}
				}
			}
		}
	}
}

void FusionFilterGive64(S** AnswerList, S* StockedE,S* OnTheFlyE,word target,int shift1,int shift2,unsigned int w,unsigned int w2,unsigned int csize){
	unsigned int i;
	S* current1;
	S* previous =NULL;
	word index = (((((*OnTheFlyE).sum[0])<<shift1)>>shift2)^target); // OnTheFlyE first draw
	if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		current1 = calloc(1,sizeof(S));
		(*current1).indice = malloc(w2*sizeof(short));
		if(SortFilter((*current1).indice,(*OnTheFlyE).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
			(*current1).sum = malloc(csize*sizeof(word));
			for (i=0; i<csize; i++){
				*((*current1).sum+i)=((*OnTheFlyE).sum[i]^(StockedE[index].sum[i]));
			}
			(*current1).next = previous;
			previous = current1;
		}
		else{	//filtered solution
			free((*current1).indice);
			free (current1);
		}
		S* current2 = &StockedE[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			current1 = calloc(1,sizeof(S));
			(*current1).indice = malloc(w2*sizeof(short));
			if(SortFilter((*current1).indice,(*OnTheFlyE).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
				(*current1).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*current1).sum+i) = ((*OnTheFlyE).sum[i]^((*current2).sum[i]));
				}
				(*current1).next = previous;
				previous = current1;
			}
			else{	//filtered solution
				free((*current1).indice);
				free (current1);
			}
		}
	}
	S* current3 = OnTheFlyE;
	while((*current3).next != NULL){ // OnTheFlyE next draws
		current3 = (*current3).next;
		word index = (((((*current3).sum[0])<<shift1)>>shift2)^target);
		if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
			current1 = calloc(1,sizeof(S));
			(*current1).indice = malloc(w2*sizeof(short));
			if(SortFilter((*current1).indice,(*current3).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
				(*current1).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*current1).sum+i)=((*current3).sum[i]^(StockedE[index].sum[i]));
				}
				(*current1).next = previous;
				previous = current1;
			}
			else{	//filtered solution
				free((*current1).indice);
				free (current1);
			}
			S* current2 = &StockedE[index];
			while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
				current2 = (*current2).next;
				current1 = calloc(1,sizeof(S));
				(*current1).indice = malloc(w2*sizeof(short));
				if(SortFilter((*current1).indice,(*current3).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
					(*current1).sum = malloc(csize*sizeof(word));
					for (i=0; i<csize; i++){
						*((*current1).sum+i)=((*current3).sum[i]^((*current2).sum[i]));
					}
					(*current1).next = previous;
					previous = current1;
				}
				else{	//filtered solution
					free((*current1).indice);
					free (current1);
				}
			}
		}
	}
	if (previous != NULL){// the solution is linked to AnswerList
		(*AnswerList)= previous;
		/*
		(*AnswerList).indice = (*previous).indice;
		(*AnswerList).sum = (*previous).sum;
		(*AnswerList).next = (*previous).next;
		*/
	}
}
#ifdef STAT
void FinalFusionFilter64(sw_list** AnswerList, S* StockedE,S* OnTheFlyE,word* Synd,int eff_word_len,unsigned int threshold,unsigned int l,unsigned int l2,unsigned int l3,unsigned int w,unsigned int w2,unsigned int csize,sw_list** observer){
#else
void FinalFusionFilter64(sw_list** AnswerList, S* StockedE,S* OnTheFlyE,word* Synd,int eff_word_len,unsigned int threshold,unsigned int l,unsigned int l2,unsigned int l3,unsigned int w,unsigned int w2,unsigned int csize){
#endif
	unsigned int i;
	int ncolumn;
	int finalweight;
	unsigned int weight;
	S builder;
	word target = ((Synd[0]<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3)));
	builder.indice = calloc(w2,sizeof(short));
	builder.sum = calloc(1,csize*sizeof(word));
	word index = ((((((*OnTheFlyE).sum[0])<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3))))^target); // OnTheFlyE first draw
	if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		if(SortFilter(builder.indice,(*OnTheFlyE).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
			/*
			printf("-");
			for (i=0;i<w2;i++){
				printf("ind %d: %d ",i,builder.indice[i]);
			}
			printf(" sum: %llx \n",((*OnTheFlyE).sum[0]^(StockedE[index].sum[0])));
			printf(" cible: %llx \n",Synd[0]);
			*/
			incr_collision_counter();
#ifdef STAT
			sw_list_add_array(observer,0,0,w2,builder.indice);
#endif
			weight = isd_weight(((*OnTheFlyE).sum[0]^(StockedE[index].sum[0]))^Synd[0]);
			if (weight < threshold) {
				ncolumn = w2;
				for (i=0;i < w2; i++){
					if (builder.indice[i] == USHRT_MAX){
						ncolumn--;
					}
				}
				incr_final_test_counter();
				final_test_cycle_stopwatch_start();
				finalweight = final_test_array(0, weight, ncolumn, builder.indice);
				final_test_cycle_stopwatch_stop();
				/*
				printf(" weight %d \n",weight);
				printf(" diff: %llx \n",((*OnTheFlyE).sum[0]^(StockedE[index].sum[0])^Synd[0]));
				*/
				if (finalweight != -1) {
					sw_list_add_array(AnswerList,0,finalweight,ncolumn,builder.indice);
				}
			}
		}
		S* current2 = &StockedE[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			if(SortFilter(builder.indice,(*OnTheFlyE).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
				/*
				printf("-");
				for (i=0;i<w2;i++){
					printf("ind %d: %d ",i,builder.indice[i]);
				}
				printf(" sum: %llx \n",((*OnTheFlyE).sum[0]^((*current2).sum[0])));
				printf(" cible: %llx \n",Synd[0]);
				*/
				incr_collision_counter();
#ifdef STAT
				sw_list_add_array(observer,0,0,w2,builder.indice);
#endif
				weight = isd_weight(((*OnTheFlyE).sum[0]^(StockedE[index].sum[0]))^Synd[0]);
				if (weight < threshold) {
					ncolumn = w2;
					for (i=0;i < w2; i++){
						if (builder.indice[i] == USHRT_MAX){
							ncolumn--;
						}
					}
					incr_final_test_counter();
					final_test_cycle_stopwatch_start();
					finalweight = final_test_array(0, weight, ncolumn, builder.indice);
					final_test_cycle_stopwatch_stop();
					/*
					printf(" weight %d \n",weight);
					printf(" diff: %llx \n",((*OnTheFlyE).sum[0]^((*current2).sum[0])^Synd[0]));
					*/
					if (finalweight != -1) {
						sw_list_add_array(AnswerList,0,finalweight,ncolumn,builder.indice);
					}
				}
			}
		}
	}
	S* current3 = OnTheFlyE;
	while((*current3).next != NULL){ // OnTheFlyE next draws
		current3 = (*current3).next;
		word index = ((((((*current3).sum[0])<<((l2+l3)+(64-eff_word_len)))>>(64 - (l-l2-l3))))^target); // OnTheFlyE first draw
		if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
			if(SortFilter(builder.indice,(*current3).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
				/*
				printf("-");
				for (i=0;i<w2;i++){
					printf("ind %d: %d ",i,builder.indice[i]);
				}
				printf(" sum: %llx \n",((*current3).sum[0]^(StockedE[index].sum[0])));
				printf(" cible: %llx \n",Synd[0]);
				*/
				incr_collision_counter();
#ifdef STAT
				sw_list_add_array(observer,0,0,w2,builder.indice);
#endif
				weight = isd_weight(((*current3).sum[0]^(StockedE[index].sum[0]))^Synd[0]);
				if (weight < threshold) {
					ncolumn = w2;
					for (i=0;i < w2; i++){
						if (builder.indice[i] == USHRT_MAX){
							ncolumn--;
						}
					}
					incr_final_test_counter();
					final_test_cycle_stopwatch_start();
					finalweight = final_test_array(0, weight, ncolumn, builder.indice);
					final_test_cycle_stopwatch_stop();
					/*
					printf(" weight %d \n",weight);
					printf(" diff: %llx \n",((*current3).sum[0]^(StockedE[index].sum[0]))^Synd[0]);
					*/
					if (finalweight != -1) {
						sw_list_add_array(AnswerList,0,finalweight,ncolumn,builder.indice);
					}
				}
			}

			S* current2 = &StockedE[index];
			while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
				current2 = (*current2).next;
				if(SortFilter(builder.indice,(*current3).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
					/*
					printf("-");
					for (i=0;i<w2;i++){
						printf("ind %d: %d ",i,builder.indice[i]);
					}
					printf(" sum: %llx \n",((*current3).sum[0]^((*current2).sum[0])));
					printf(" cible: %llx \n",Synd[0]);
					*/
					incr_collision_counter();
#ifdef STAT
					sw_list_add_array(observer,0,0,w2,builder.indice);
#endif
					weight = isd_weight(((*current3).sum[0]^((*current2).sum[0]))^Synd[0]);
					if (weight < threshold) {
						ncolumn = w2;
						for (i=0;i < w2; i++){
							if (builder.indice[i] == USHRT_MAX){
								ncolumn--;
							}
						}
						incr_final_test_counter();
						final_test_cycle_stopwatch_start();
						finalweight = final_test_array(0, weight, ncolumn, builder.indice);
						final_test_cycle_stopwatch_stop();
						/*
						printf(" weight %d \n",weight);
						printf(" diff: %llx \n",((*current3).sum[0]^((*current2).sum[0]))^Synd[0]);
						*/
						if (finalweight != -1) {
							sw_list_add_array(AnswerList,0,finalweight,ncolumn,builder.indice);
						}
					}
				}
			}
		}
	}
	free(builder.sum);
	free(builder.indice);
}

void Sort(unsigned short* dest,unsigned short* s1,unsigned short* s2,unsigned short size1,unsigned short size2){
	unsigned short p1 =0;
	unsigned short p2 =0;
	while (p1 < size1 && p2 < size2){ //sorting
		if (*(s1+p1) < *(s2+p2) ){
			*(dest+p1+p2)= *(s1+p1);
			p1++;
		}
		else{
			*(dest+p1+p2)= *(s2+p2);
			p2++;
		}
	}
	while(p1<size1){ // filling the end of dest. Only one of the two while loop will do something
		*(dest+p1+p2)= *(s1+p1);
		p1++;
	}
	while(p2<size2){ // filling the end of dest. Only one of the two while loop will do something
		*(dest+p1+p2)= *(s2+p2);
		p2++;
	}
}

int SortFilter(unsigned short* dest,unsigned short* s1,unsigned short* s2,unsigned short size1,unsigned short size2,unsigned short targetsize){
	unsigned short p1 =0;
	unsigned short p2 =0;
	unsigned short currentsize = 1;
	int i;
	if (*(s1+p1) < *(s2+p2)){ // first step
		*dest= *(s1+p1);
		p1++;
	}
	else{
		*dest= *(s2+p2);
		p2++;
	}
	while (p1 < size1 && p2 < size2){ //sorting
		if (*(s1+p1) < *(s2+p2) ){
			if (*(dest+currentsize-1) == *(s1+p1)){
				p1++;
			}
			else {
				currentsize++;
				if(currentsize > targetsize){
					return 0; //filtered solution
				}
				*(dest+currentsize-1)= *(s1+p1);
				p1++;
			}
		}
		else{
			if ((*(dest+currentsize-1) == *(s2+p2)) || (*(s2+p2) == USHRT_MAX)){
				p2++;
			}
			else {

				currentsize++;
				if(currentsize > targetsize){
					return 0; //filtered solution
				}
				*(dest+currentsize-1)= *(s2+p2);
				p2++;
			}
		}
	}
	while(p1<size1){ // filling the end of dest. Only one of the two while loop will do something
		if ((*(dest+currentsize-1) == *(s1+p1)) || (*(s1+p1) == USHRT_MAX)){
			p1++;
		}
		else {
			currentsize++;
			if(currentsize > targetsize){
				return 0; //filtered solution
			}
			*(dest+currentsize-1)= *(s1+p1);
			p1++;
		}
	}
	while(p2<size2){ // filling the end of dest. Only one of the two while loop will do something
		if ((*(dest+currentsize-1) == *(s2+p2)) || (*(s2+p2) == USHRT_MAX)){
			p2++;
		}
		else {
			currentsize++;
			if(currentsize > targetsize){
				return 0; //filtered solution
			}
			*(dest+currentsize-1)= *(s2+p2);
			p2++;
		}
	}
	for (i=currentsize; i<targetsize ; i++){ //potential garbage erasing
		*(dest+i) = USHRT_MAX;
	}
	return 1;
}

void freelist(S* draw){
	if((*draw).next != NULL){
		freelist((*draw).next);
		free((*draw).next);
		(*draw).next= NULL;
	}
	free((*draw).sum);
	(*draw).sum= NULL;
	free((*draw).indice);
	(*draw).indice= NULL;
}

void PrintDoubleStat(sw_list** AnswerList,unsigned short p){

	int sorter(const void* a, const void* b){
		unsigned int i;
		for (i=0;i<p;i++){
			if ((*((sw_list*) a)).pos[i] < (*((sw_list*) b)).pos[i]){
				return -1;
			}
			else if ((*((sw_list*) a)).pos[i] > (*((sw_list*) b)).pos[i]){
				return 1;
			}
		}
		return 0;
	}

	int stater(sw_list* list,unsigned int i,unsigned int j,unsigned int* stat,unsigned int p,unsigned int size){
		unsigned int l;
		if(list[i].pos[p-1]==list[i-1].pos[p-1]){
			if (i<size-1){
				l=stater(list,i+1,j+1,stat,p,size);
				if (l==0){
					if (j>19){
						stat[19]++;
					}
					else{
						stat[j]++;
					}
				}
				return l+1;
			}
			else {
				return 0;
			}
		}
		else {
			return 0;
		}
	}

	unsigned int i;
	unsigned int stat[20];
	unsigned int size=0;
	int inc;
	double prop[20];
	for (i=0;i<20;i++){
		prop[i]=0.0;
		stat[i]=0;
	}
	sw_list* current = *AnswerList;
	sw_list* list;
	while (current !=NULL){
		size++;
		current=(*current).next;
	}
	list=calloc(size,sizeof(sw_list));
	current = *AnswerList;
	for (i=0;i<size;i++){
		//if ((*current).weight==p){
		list[i]=*current;
		//}
		current=(*current).next;

	}
	qsort(list,size,sizeof(sw_list),sorter);

	i=1;
	while(i<size-2){
		inc = stater(list,i,0,stat,p,size);
		i+=inc+1;
	}
	printf("iteration a %d collisions) \n",size);
	for (i=0;i<19;i++){
	prop[i] = ((double) stat[i])*(i+1) / ((double) size) * 100.0;
	printf("pourcentage de %d-uplet: %f %% \n",i+2,prop[i]);
	}
	prop[19] = ((double) stat[19])*(21) / ((double) size) * 100.0;
	printf("pourcentage de N-uplet superieur a %d: %f %% environ (minorantt) \n",21,prop[19]);

	free(list);
}







/*
 * debug version
 *
 * Use to see EStep1 elements
 *
void fusionstore1(S* EStep1,word target,S* table,int shift1,int shift2,int shift3,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	word sumr1;
	word index = ((sums[w/2-1]>>shift1)^target);
	printf("new store \n");
	fflush(stdout);
	if (table[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		sumr1 = (((sums[w/2-1]^(table[index].sum[0]))<<shift2)>>shift3);
		if(EStep1[sumr1].indice==NULL){			//there is not yet a solution at sumr1
			printf("first sol");
			fflush(stdout);
			EStep1[sumr1].sum = malloc(csize*sizeof(word));
			EStep1[sumr1].indice = malloc(w*sizeof(short));
			for (i=0; i<csize; i++){
				*(EStep1[sumr1].sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
			}
			Sort(EStep1[sumr1].indice,table[index].indice,indice,w/2,w/2);
			printf(" mempos: %u  \n",(&EStep1[sumr1]));
			fflush(stdout);
			printf(" next: %u  \n",EStep1[sumr1].next);
			fflush(stdout);
			printf("indice: %d \n",(EStep1[sumr1].indice[0]));
			fflush(stdout);
			printf("indice: %d \n",(EStep1[sumr1].indice[1]));
			fflush(stdout);
			printf("sum: %u \n",(unsigned int) ((EStep1[sumr1].sum[0])>>32));
			fflush(stdout);
			gets(stdin);
		}
		else{			//there are already one or more solutions at sumr1
			S* current1 = &EStep1[sumr1];
			printf("other sol");
			fflush(stdout);
			while((*current1).next!=NULL){
				current1 = (*current1).next;
			}
			(*current1).next= calloc(1,sizeof(S));
			(*((*current1).next)).indice = malloc(w*sizeof(short));
			Sort((*((*current1).next)).indice,table[index].indice,indice,w/2,w/2);
			(*((*current1).next)).sum = malloc(csize*sizeof(word));
			for (i=0; i<csize; i++){
				*((*((*current1).next)).sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
			}
		printf(" mempos: %u  \n",((*current1).next));
		fflush(stdout);
		printf(" next: %u  \n",(*((*current1).next)).next);
		fflush(stdout);
		printf("indice: %d \n",((*((*current1).next)).indice[0]));
		fflush(stdout);
		printf("indice: %d \n",((*((*current1).next)).indice[1]));
		fflush(stdout);
		printf("sum: %u \n",(unsigned int)(((*((*current1).next)).sum[0])>>32));
		fflush(stdout);
		gets(stdin);
		}

		S* current2 = &table[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			sumr1 = (((sums[w/2-1]^((*current2).sum[0]))<<shift2)>>shift3);

			if(EStep1[sumr1].indice==NULL){			//there is not yet a solution at sumr1
				printf("first sol");
				fflush(stdout);
				EStep1[sumr1].sum = malloc(csize*sizeof(word));
				EStep1[sumr1].indice = malloc(w*sizeof(short));
				for (i=0; i<csize; i++){
					*(EStep1[sumr1].sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
				}
				Sort(EStep1[sumr1].indice,(*current2).indice,indice,w/2,w/2);
				printf(" mempos: %u  \n",(&EStep1[sumr1]));
				fflush(stdout);
				printf(" next: %u  \n",EStep1[sumr1].next);
				fflush(stdout);
				printf("indice: %d \n",(EStep1[sumr1].indice[0]));
				fflush(stdout);
				printf("indice: %d \n",(EStep1[sumr1].indice[1]));
				fflush(stdout);
				printf("sum: %u \n",(unsigned int) ((EStep1[sumr1].sum[0])>>32));
				fflush(stdout);
				gets(stdin);
			}
			else{			//there are already one or more solutions at sumr1
				printf("other sol");
				fflush(stdout);
				S* current1 = &EStep1[sumr1];
				while((*current1).next!=NULL){
					current1 = (*current1).next;
				}
				(*current1).next= calloc(1,sizeof(S));
				(*((*current1).next)).indice = malloc(w*sizeof(short));
				Sort((*((*current1).next)).indice,(*current2).indice,indice,w/2,w/2);
				(*((*current1).next)).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*((*current1).next)).sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
				}
			printf(" mempos: %u  \n",((*current1).next));
			fflush(stdout);
			printf(" next: %u  \n",(*((*current1).next)).next);
			fflush(stdout);
			printf("indice: %d \n",((*((*current1).next)).indice[0]));
			fflush(stdout);
			printf("indice: %d \n",((*((*current1).next)).indice[1]));
			fflush(stdout);
			printf("sum: %u \n",(unsigned int) (((*((*current1).next)).sum[0])>>32));
			fflush(stdout);
			gets(stdin);
			}

		}
	}
}
*/

/*
 * debug version
 *
 * Use to see given elements
 *
void fusiongive1(S* answer,word target,S* table,int shift1,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
	unsigned int i;
	printf(" new give \n");
	fflush(stdout);
	word index = ((sums[w/2-1]>>shift1)^target);
	if (table[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		(*answer).sum = malloc(csize*sizeof(word));
		(*answer).indice = malloc(w*sizeof(short));
		for (i=0; i<csize; i++){
			*((*answer).sum+i)=(sums[(w/2)*(1+i)-1]^(table[index].sum[i]));
		}
		Sort((*answer).indice,table[index].indice,indice,w/2,w/2);
		printf(" mempos: %u  \n",answer);
		fflush(stdout);
		printf(" next: %u  \n",((*answer).next));
		fflush(stdout);
		printf("indice: %d \n",((*answer).indice[0]));
		fflush(stdout);
		printf("indice: %d \n",((*answer).indice[1]));
		fflush(stdout);
		//printf("indice: %d \n",((*answer).indice[2]));
		//fflush(stdout);
		//printf("indice: %d \n",((*answer).indice[3]));
		//fflush(stdout);
		printf("sum: %u \n",(unsigned int)(((*answer).sum[0])>>32));
		fflush(stdout);
		gets(stdin);
		S* current1 = answer;
		S* current2 = &table[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			(*current1).next=calloc(1,sizeof(S));
			current1 = (*current1).next;
			(*current1).sum = malloc(csize*sizeof(word));
			(*current1).indice = malloc(w*sizeof(short));
			for (i=0; i<csize; i++){
				*((*current1).sum+i)=(sums[(w/2)*(1+i)-1]^((*current2).sum[i]));
			}
			Sort((*current1).indice,(*current2).indice,indice,w/2,w/2);
			printf(" mempos: %u  \n",current1);
			fflush(stdout);
			printf(" next: %u  \n",((*current1).next));
			fflush(stdout);
			printf("indice: %d \n",((*current1).indice[0]));
			fflush(stdout);
			printf("indice: %d \n",((*current1).indice[1]));
			fflush(stdout);
			//printf("indice: %d \n",((*current1).indice[2]));
			//fflush(stdout);
			//printf("indice: %d \n",((*current1).indice[3]));
			//fflush(stdout);
			printf("sum: %u \n",(unsigned int)(((*current1).sum[0])>>32));
			fflush(stdout);
			gets(stdin);
		}
	}
}
*
*/
