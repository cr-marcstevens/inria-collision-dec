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

void fusionstore1(S* EStep1,word target,S* table,int shift1,int shift2,int shift3,unsigned short* indice,word* sums,unsigned int w,unsigned int csize){
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
			sumr1 = (((sums[w-1]^((*current2).sum[0]))<<shift2)>>shift3);

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
	word index = ((sums[w-1]>>shift1)^target);
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
	word index = (((((*OnTheFlyE).sum[0])^target)<<shift1)>>shift2); // OnTheFlyE first draw
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


		S* current2 = &StockedE[index];
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
		word index = (((((*current3).sum[0])^target)<<shift1)>>shift2);
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

void FusionFilterGive64(S* AnswerList, S* StockedE,S* OnTheFlyE,word target,int shift1,int shift2,unsigned int w,unsigned int w2,unsigned int csize){
	unsigned int i;
	S* current1 = AnswerList;
	word index = (((((*OnTheFlyE).sum[0])^target)<<shift1)>>shift2); // OnTheFlyE first draw
	if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
		(*current1).indice = malloc(w2*sizeof(short));
		if(SortFilter((*current1).indice,(*OnTheFlyE).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
			(*current1).sum = malloc(csize*sizeof(word));
			for (i=0; i<csize; i++){
				*((*current1).sum+i)=((*OnTheFlyE).sum[i]^(StockedE[index].sum[i]));
			}
			(*current1).next = calloc(1,sizeof(S));
			current1 = (*current1).next;
		}
		else{	//filtered solution
			free((*current1).indice);
			(*current1).indice = NULL;
		}
		S* current2 = &StockedE[index];
		while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
			current2 = (*current2).next;
			(*current1).indice = malloc(w2*sizeof(short));
			if(SortFilter((*current1).indice,(*OnTheFlyE).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
				(*current1).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*current1).sum+i) = ((*OnTheFlyE).sum[i]^((*current2).sum[i]));
				}
				(*current1).next = calloc(1,sizeof(S));
				current1 = (*current1).next;
			}
			else{	//filtered solution
				free((*current1).indice);
				(*current1).indice = NULL;
			}
		}
	}
	S* current3 = OnTheFlyE;
	while((*current3).next != NULL){ // OnTheFlyE next draws
		current3 = (*current3).next;
		word index = (((((*current3).sum[0])^target)<<shift1)>>shift2);
		if (StockedE[index].indice != NULL){		// there is a corresponding solution to build the targeted syndrome
			(*current1).indice = malloc(w2*sizeof(short));
			if(SortFilter((*current1).indice,(*current3).indice,StockedE[index].indice,w,w,w2)){ //The solution pass the filter
				(*current1).sum = malloc(csize*sizeof(word));
				for (i=0; i<csize; i++){
					*((*current1).sum+i)=((*current3).sum[i]^(StockedE[index].sum[i]));
				}
				(*current1).next = calloc(1,sizeof(S));
				current1 = (*current1).next;
			}
			else{	//filtered solution
				free((*current1).indice);
				(*current1).indice = NULL;
			}
			S* current2 = &StockedE[index];
			while((*current2).next != NULL){		// there are other corresponding solutions to build the targeted syndrome
				current2 = (*current2).next;
				(*current1).indice = malloc(w2*sizeof(short));
				if(SortFilter((*current1).indice,(*current3).indice,(*current2).indice,w,w,w2)){ //The solution pass the filter
					(*current1).sum = malloc(csize*sizeof(word));
					for (i=0; i<csize; i++){
						*((*current1).sum+i)=((*current3).sum[i]^((*current2).sum[i]));
					}
					(*current1).next = calloc(1,sizeof(S));
					current1 = (*current1).next;
				}
				else{	//filtered solution
					free((*current1).indice);
					(*current1).indice = NULL;
				}
			}
		}
	}
}

void FinalFusionFilter64(){

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
				*(dest+currentsize)= *(s1+p1);
				p1++;
				currentsize++;
				if(currentsize > targetsize){
					return 0; //filtered solution
				}
			}
		}
		else{
			if (*(dest+currentsize-1) == *(s2+p2)){
				p2++;
			}
			else {
				*(dest+currentsize)= *(s2+p2);
				p2++;
				currentsize++;
				if(currentsize > targetsize){
					return 0; //filtered solution
				}
			}
		}
	}
	while(p1<size1){ // filling the end of dest. Only one of the two while loop will do something
		if (*(dest+currentsize-1) == *(s1+p1)){
			p1++;
		}
		else {
			*(dest+currentsize)= *(s1+p1);
			p1++;
			currentsize++;
			if(currentsize > targetsize){
				return 0; //filtered solution
			}
		}
	}
	while(p2<size2){ // filling the end of dest. Only one of the two while loop will do something
		if (*(dest+currentsize-1) == *(s2+p2)){
			p2++;
		}
		else {
			*(dest+currentsize)= *(s2+p2);
			p2++;
			currentsize++;
			if(currentsize > targetsize){
				return 0; //filtered solution
			}
		}
	}
	return 1;
}

void freelist(S draw){
	if(draw.next != NULL){
		freelist(*(draw.next));
		free(draw.next);
	}
	free(draw.sum);
	free(draw.indice);
}

