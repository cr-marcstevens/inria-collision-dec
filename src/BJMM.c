/*
 * BJMM.c
 *
 *  Created on: 25 juin 2013
 *      Author: Mathieu ARIA
 */

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "sub_isd.h"
#include <m4ri/m4ri.h>
#include "libisd.h"
#include "sparse_words_list.h"
#include "final_test.h"
#include "measure.h"
#include "cpucycles/cpucycles.h"
#include "BJMMtools.h"
#include "measure.h"
#include "prng.h"

static ranctx* state;
static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, l2, l3, p, e1, e2, w, L_len, threshold, csize;
static unsigned int p2;
static unsigned int p1;
static word s1, s2, s3, s4, a, sma;
static int eff_word_len,shift1,shift2,shift1p;
static sw_list** h;
#ifdef STAT
static sw_list** observer;
#endif

void print_parameters(isd_params* params) {
	printf("n : %d\n", params->n);
	printf("r : %d\n", params->r);
	printf("w : %d\n", params->w);
	printf("l : %d\n", params->l);
	printf("l2 : %d\n", params->l2);
	printf("l3 : %d\n", params->l3);
	printf("e1 : %d\n", params->e1);
	printf("e2 : %d\n", params->e2);
	printf("p : %d\n", params->p);
	printf("eff_word_len : %ld\n", min(params->r, word_len));
	printf("threshold : %d\n", params->weight_threshold);
}

void sub_isd_init(isd_params* params, word* local_L, word* local_synds, unsigned int local_N, sw_list** local_h, ranctx* local_state) {
	//printf("initiate BJMM \n");
	//fflush(stdout);
	print_parameters(params);

	L = local_L;
	N = local_N;
	syndsprime = local_synds;
	n = params->n;
	r = params->r;
	k = params->k;
	l = params->l;
	l2 = params->l2; //  named r2 in the paper
	l3 = params->l3; //  named r1 in the paper
	e1 = params->e1; //  filter 1 parameter
	e2 = params->e2; //  filter 2 parameter
	p = params->p;
	w = params->w;
	h = local_h;
	state = local_state;
	csize = params->csize;


	L_len = k+l;

	threshold = params->weight_threshold;
	eff_word_len = min(r, word_len);
	shift1 = eff_word_len - l2;
	shift1p = l2 + (64 - eff_word_len);
	shift2 = 64 - l3;

	p1= p/2 + e1;
	p2 = p1/2 + e2;

	if((p2/2)!=((p2+1)/2)){
		printf("WARNING: odd p2 not handled \n");
		fflush(stdout);
	}

	// hard-coded parameters. Imply l2 >= 2
	//s1 = 0;
	//s2 = 1;
	//s3 = 2;
	//a = 1;

	//printf("BJMM initiated \n");
	//fflush(stdout);
}


void sub_isd() {
#ifdef STAT
	observer=calloc(1,sizeof(sw_list*));
	*observer=NULL;
#endif

	s1 = random_range(state, 1UL<<l2);
	s2 = random_range(state, 1UL<<l2);
	s3 = random_range(state, 1UL<<l2);
	s4 = (s1^s2^s3)^(syndsprime[0] >> shift1);

	//printf("synd: %u \n",(unsigned int)((syndsprime[0])>>shift1));
	//printf("s1234: %u \n",(unsigned int)((s1^s2^s3^s4)));


	sma = (((syndsprime[0])<<shift1p)>>shift2)^a;
	/*
	 * Lists construction from E1, E2, E3, E4
	 *
	 * warning: in the whole function, we assume l2+l3 <= 64)
	 */
	unsigned int i,ii,j;
	unsigned short* indice = malloc((p2/2)*sizeof(unsigned short));
	word* sums = malloc((p2/2)*csize*sizeof(word));
	S* lists[4];
	S* EStep1[2];
	S* EStep2;
	S* temp;
	S** temp2 = calloc(1,sizeof(S*));

	temp = calloc(1,sizeof(S));
	//temp2 = calloc(1,sizeof(S));
	for (i=0; i<4;i++){
		//printf("building list %i \n",i);
		lists[i] = calloc((1ULL << l2),sizeof(S));
		/* initialization of indice
		 *  We assume p2/2 << k+l.
		 *  Then the p2/2 first colums will always be in the first block in all distributions
		 *  (except the second which is odd/even).
		 *  Please choose (k+l) divisible by 4 to avoid possible bugs
		 *  due to blocks of non equal size in distribution 3 and 4.
		 */
		switch (i){
			case 0:
				for (j=0;j<(p2/2);j++){
					indice[j]=j;
				} break;
			case 1:
				for (j=0;j<(p2/2);j++){
					indice[j]=2*j;
				} break;

			case 2:
				for (j=0;j<(p2/2);j++){
					indice[j]=j;
				} break;
			case 3:
				for (j=0;j<(p2/2);j++){
					indice[j]=j;
				} break;
		}
		// initialization of sums
		for (ii=0;ii<csize;ii++){
			sums[ii*(p2/2)]=L[indice[0]+ii*L_len];
			for (j=1;j<(p2/2);j++){
				sums[j+ii*(p2/2)]=sums[(j-1)+ii*(p2/2)]^L[(indice[j]+ii*L_len)];
			}
		}
		//building loop
		h1store(lists[i],shift1,indice,sums,(p2/2),csize);
		while(next(L,indice,sums,csize,L_len,(p2/2),i)){
				h1store(lists[i],shift1,indice,sums,(p2/2),csize);
		}
	}

	// Building all Ei lists on the fly.
	//---------------------

	//printf("fusion -> E1 \n");
	EStep1[0] = calloc((1ULL << l3),sizeof(S));
	for (j=0;j<(p2/2);j++){
		indice[j]=(L_len/2)+j;
	}
	for (ii=0;ii<csize;ii++){
		sums[ii*(p2/2)]=L[indice[0]+ii*L_len];
		for (j=1;j<(p2/2);j++){
			sums[j+ii*(p2/2)]=sums[(j-1)+ii*(p2/2)]^L[indice[j]+ii*L_len];
		}
	}
	fusionstore1(EStep1[0],s1,lists[0],shift1,shift1p,shift2,indice,sums,p2,csize);
	while(next2(L,indice,sums,csize,L_len,(p2/2),0)){
		fusionstore1(EStep1[0],s1,lists[0],shift1,shift1p,shift2,indice,sums,p2,csize);
	}
	//used list deletion
	for (j=0;j<(1UL << l2);j++){
		if (lists[0][j].indice != NULL){
			freelist(&(lists[0][j]));
		}
	}
	free(lists[0]);
	//---------------------
	if ((l)<=64) {
		//printf("fusion ->E2 and E1 + E2 -> E'1 \n");
		for (j=0;j<(p2/2);j++){
		indice[j]=(2*j)+1;
		}
		for (ii=0;ii<csize;ii++){
			sums[ii*(p2/2)]=L[indice[0]+ii*L_len];
			for (j=1;j<(p2/2);j++){
				sums[j+ii*(p2/2)]=sums[(j-1)+ii*(p2/2)]^L[indice[j]+ii*L_len];
			}
		}
		EStep2 = calloc((1ULL << (l-l2-l3)),sizeof(S));

		fusiongive1(temp,s2,lists[1],shift1,indice,sums,p2,csize);
		if ((*temp).indice != NULL){
			FusionFilterStore64(EStep2,EStep1[0],temp,a,shift1p,shift2,eff_word_len,l,l2,l3,p2,p1,csize);
			freelist(temp);
			(*temp).indice = NULL;
		}
			while(next2(L,indice,sums,csize,L_len,(p2/2),1)){
				fusiongive1(temp,s2,lists[1],shift1,indice,sums,p2,csize);
				if ((*temp).indice != NULL){
					//printf("before ");
					//fflush(stdout);
					FusionFilterStore64(EStep2,EStep1[0],temp,a,shift1p,shift2,eff_word_len,l,l2,l3,p2,p1,csize);
					//printf(" after ");
					//fflush(stdout);
					freelist(temp);
					//printf("ok \n");
					//fflush(stdout);
					(*temp).indice = NULL;
				}
			}


		//used list deletion
		for (j=0;j<(1UL << l2);j++){
			if (lists[1][j].indice != NULL){
				freelist(&(lists[1][j]));
			}
		}
		free(lists[1]);
		//used list deletion
		for (j=0;j<(1UL << l3);j++){
			if (EStep1[0][j].indice != NULL){
				freelist(&(EStep1[0][j]));
			}
		}
		free(EStep1[0]);
		//---------------------
		//printf("fusion -> E3 \n");
		EStep1[1]= calloc((1UL << l3),sizeof(S));
		for (j=0;j<(p2/2);j++){
			indice[j]=(L_len/4)+j;
		}
		for (ii=0;ii<csize;ii++){
			sums[ii*(p2/2)]=L[indice[0]+ii*L_len];
			for (j=1;j<(p2/2);j++){
				sums[j+ii*(p2/2)]=sums[(j-1)+ii*(p2/2)]^L[indice[j]+ii*L_len];
			}
		}
		fusionstore1(EStep1[1],s3,lists[2],shift1,shift1p,shift2,indice,sums,p2,csize);
		while(next2(L,indice,sums,csize,L_len,(p2/2),2)){
			fusionstore1(EStep1[1],s3,lists[2],shift1,shift1p,shift2,indice,sums,p2,csize);
		}
		//used list deletion
		for (j=0;j<(1UL << l2);j++){
			if (lists[2][j].indice != NULL){
				freelist(&(lists[2][j]));
			}
		}
		free(lists[2]);
		//---------------------
		//printf("fusion ->E4 and E3 + E4 -> E'2 and E'1 + E'2 -> E \n");
		for (j=0;j<(p2/2);j++){
			indice[j]=(L_len/4)+j;
		}
		for (ii=0;ii<csize;ii++){
			sums[ii*(p2/2)]=L[indice[0]+ii*L_len];
			for (j=1;j<(p2/2);j++){
				sums[j+ii*(p2/2)]=sums[(j-1)+ii*(p2/2)]^L[indice[j]+ii*L_len];
			}
		}
		fusiongive1(temp,s4,lists[3],shift1,indice,sums,p2,csize);
		if ((*temp).indice != NULL){
			FusionFilterGive64(temp2,EStep1[1],temp,sma,shift1p,shift2,p2,p1,csize);
			if ( (*temp2) != NULL){
#ifdef STAT
				FinalFusionFilter64(h,EStep2,*temp2,syndsprime,eff_word_len,threshold,l,l2,l3,p1,p,csize,observer);
#else
				FinalFusionFilter64(h,EStep2,*temp2,syndsprime,eff_word_len,threshold,l,l2,l3,p1,p,csize);
#endif
				freelist(*temp2);
				//(*temp2).indice = NULL;
				free(*temp2);
				*temp2= NULL;
			}
			freelist(temp);
			(*temp).indice = NULL;
		}



		while(next2(L,indice,sums,csize,L_len,(p2/2),3)){
			fusiongive1(temp,s4,lists[3],shift1,indice,sums,p2,csize);
			if ((*temp).indice != NULL){
				FusionFilterGive64(temp2,EStep1[1],temp,sma,shift1p,shift2,p2,p1,csize);
				if ( (*temp2) != NULL){
#ifdef STAT
					FinalFusionFilter64(h,EStep2,*temp2,syndsprime,eff_word_len,threshold,l,l2,l3,p1,p,csize,observer);
#else
					FinalFusionFilter64(h,EStep2,*temp2,syndsprime,eff_word_len,threshold,l,l2,l3,p1,p,csize);
#endif
					freelist(*temp2);
					//(*temp2).indice = NULL;
					free(*temp2);
					*temp2= NULL;
				}
				freelist(temp);
				(*temp).indice = NULL;
			}
		}

		//used list deletion
		for (j=0;j<(1UL << l2);j++){
			if (lists[3][j].indice != NULL){
				freelist(&(lists[3][j]));
			}
		}
		free(lists[3]);
		//used list deletion
		for (j=0;j<(1UL << l3);j++){
			if (EStep1[1][j].indice != NULL){
				freelist(&(EStep1[1][j]));
			}
		}
		free(EStep1[1]);
		//used list deletion
		for (j=0;j<(1UL << (l-l2-l3));j++){
			if (EStep2[j].indice != NULL){
				freelist(&(EStep2[j]));
			}
		}
		free(EStep2);
	}
	else {
		//TODO l > 64 not implemented!
		printf("warning! l > 64 not implemented \n");
	}
	free(temp);
	free(temp2);
	free(indice);
	free(sums);
	// BJMM done
	//printf("-");
	//printf("h: %llu", h[0]);

#ifdef STAT
	PrintDoubleStat(observer,p);
	sw_list_free(*observer);
	free(observer);
#endif

}

void sub_isd_free() {

}
