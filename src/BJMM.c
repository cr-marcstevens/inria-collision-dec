/*
 * BJMM.c
 *
 *  Created on: 25 juin 2013
 *      Author: Mathieu ARIA
 */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "sub_isd.h"
#include "m4ri/m4ri.h"
#include "libisd.h"
#include "sparse_words_list.h"
#include "final_test.h"
#include "measure.h"
#include "cpucycles/cpucycles.h"
#include "BJMMtools.h"


static word* L;
static unsigned int N;
static word* syndsprime;
static unsigned int n, k, r, l, l2, l3, p, e1, e2, w, L_len, threshold, csize;
static unsigned int p2;
static unsigned int p1;
static word s1, s2, s3, s4, a, sma;
static int eff_word_len,shift1,shift2,shift1p;
static sw_list** h;

void sub_isd_init(word* simple_HprimemodT, unsigned int local_N, word* local_syndsprime, unsigned int local_n, unsigned int local_r,unsigned int local_l, unsigned int local_l2, unsigned int local_l3, unsigned int local_p, unsigned int local_e1, unsigned int local_e2, unsigned int local_w, unsigned int local_threshold,unsigned int local_csize, sw_list** local_h) {
	L = simple_HprimemodT;
	N = local_N;
	syndsprime = local_syndsprime;
	n = local_n;
	r = local_r;
	l = local_l;
	l2 = local_l2; //  named r2 in the paper
	l3 = local_l3; //  named r1 in the paper
	e1 = local_e1; //  filter 1 parameter
	e2 = local_e2; //  filter 2 parameter
	p = local_p;
	w = local_w;
	h = local_h;
	csize= local_csize;

	k = n-r;

	L_len = k+l;

	threshold = local_threshold;
	eff_word_len = min(r, word_len);
	shift1 = eff_word_len - l2;
	shift1p = l2 + (64 - eff_word_len);
	shift2 = 64 - l3;

	p1= p/2 + e1;
	p2 = p1/2 + e2;

	// hard-coded parameters. Imply l2 >= 2
	s1 = 0;
	s2 = 1;
	s3 = 2;
	a = 1;

}


void sub_isd() {

	s4 = (s1^s2^s3)^(syndsprime[0] >> shift1);
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
	S* temp2;

	temp = calloc(1,sizeof(S));
	temp2 = calloc(1,sizeof(S));
	for (i=0; i<4;i++){
		printf("building list %i",i);
		lists[i] = calloc((1UL << l2),sizeof(S));
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
			sums[ii*w]=L[indice[0]+ii*L_len];
			for (j=1;j<w;j++){
				sums[j+ii*w]=sums[(j-1)+ii*w]^L[indice[j]+ii*L_len];
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
	printf("fusion -> E1");
	EStep1[1] = calloc((1UL << l3),sizeof(S));
	for (j=0;j<(p2/2);j++){
		indice[j]=(L_len/2)+j;
	}
	for (ii=0;ii<csize;ii++){
		sums[ii*w]=L[indice[0]+ii*L_len];
		for (j=1;j<w;j++){
			sums[j+ii*w]=sums[(j-1)+ii*w]^L[indice[j]+ii*L_len];
		}
	}
	fusionstore1(EStep1[1],s1,lists[1],shift1,shift1p,shift2,indice,sums,p2,csize);
	while(next2(L,indice,sums,csize,L_len,(p2/2),1)){
		fusionstore1(EStep1[1],s1,lists[1],shift1,shift1p,shift2,indice,sums,p2,csize);
	}
	//used list deletion
	for (j=0;j<(1UL << l2);j++){
		if (lists[0][j].indice != NULL){
			freelist(lists[0][j]);
		}
	}
	free(lists[0]);
	//---------------------
	if ((l)<=64) {
		printf("fusion ->E2 and E1 + E2 -> E'1");
		for (j=0;j<(p2/2);j++){
		indice[j]=(2*j)+1;
		}
		for (ii=0;ii<csize;ii++){
			sums[ii*w]=L[indice[0]+ii*L_len];
			for (j=1;j<w;j++){
				sums[j+ii*w]=sums[(j-1)+ii*w]^L[indice[j]+ii*L_len];
			}
		}
		EStep2 = calloc((1UL << (l-l2-l3)),sizeof(S));

		fusiongive1(temp,s2,lists[2],shift1,indice,sums,p2,csize);
		FusionFilterStore64(EStep2,EStep1[1],temp,a,shift1p,shift2,eff_word_len,l,l2,l3,p2,p1,csize);
		freelist(*temp);
		while(next2(L,indice,sums,csize,L_len,(p2/2),2)){
			fusiongive1(temp,s2,lists[2],shift1,indice,sums,p2,csize);
			FusionFilterStore64(EStep2,EStep1[1],temp,a,shift1p,shift2,eff_word_len,l,l2,l3,p2,p1,csize);
			freelist(*temp);
		}

		//used list deletion
		for (j=0;j<(1UL << l2);j++){
			if (lists[1][j].indice != NULL){
				freelist(lists[1][j]);
			}
		}
		free(lists[1]);
		//used list deletion
		for (j=0;j<(1UL << l3);j++){
			if (EStep1[0][j].indice != NULL){
				freelist(EStep1[0][j]);
			}
		}
		free(EStep1[0]);
		//---------------------
		printf("fusion -> E3");
		EStep1[2]= calloc((1UL << l3),sizeof(S));
		for (j=0;j<(p2/2);j++){
			indice[j]=(L_len/4)+j;
		}
		for (ii=0;ii<csize;ii++){
			sums[ii*w]=L[indice[0]+ii*L_len];
			for (j=1;j<w;j++){
				sums[j+ii*w]=sums[(j-1)+ii*w]^L[indice[j]+ii*L_len];
			}
		}
		fusionstore1(EStep1[2],s3,lists[3],shift1,shift1p,shift2,indice,sums,p2,csize);
		while(next2(L,indice,sums,csize,L_len,(p2/2),3)){
			fusionstore1(EStep1[2],s3,lists[3],shift1,shift1p,shift2,indice,sums,p2,csize);
		}
		//used list deletion
		for (j=0;j<(1UL << l2);j++){
			if (lists[2][j].indice != NULL){
				freelist(lists[2][j]);
			}
		}
		free(lists[2]);
		//---------------------
		printf("fusion ->E4 and E3 + E4 -> E'2 and E'1 + E'2 -> E ");
		for (j=0;j<(p2/2);j++){
			indice[j]=(L_len/4)+j;
		}
		for (ii=0;ii<csize;ii++){
			sums[ii*w]=L[indice[0]+ii*L_len];
			for (j=1;j<w;j++){
				sums[j+ii*w]=sums[(j-1)+ii*w]^L[indice[j]+ii*L_len];
			}
		}
		fusiongive1(temp,s2,lists[2],shift1,indice,sums,p2,csize);
		FusionFilterGive64(temp2,EStep1[2],temp,sma,shift1p,shift2,p2,p1,csize);
		FinalFusionFilter64(h[0],EStep2,temp2,syndsprime,eff_word_len,l,l2,l3,p1,p,csize);
		freelist(*temp);
		freelist(*temp2);
		while(next2(L,indice,sums,csize,L_len,(p2/2),2)){
			fusiongive1(temp,s2,lists[2],shift1,indice,sums,p2,csize);
			FusionFilterGive64(temp2,EStep1[2],temp,sma,shift1p,shift2,p2,p1,csize);
			FinalFusionFilter64(h[0],EStep2,temp2,syndsprime,eff_word_len,l,l2,l3,p,p,csize);
			freelist(*temp);
			freelist(*temp2);
		}

		//used list deletion
		for (j=0;j<(1UL << l2);j++){
			if (lists[3][j].indice != NULL){
				freelist(lists[3][j]);
			}
		}
		free(lists[3]);
		//used list deletion
		for (j=0;j<(1UL << l3);j++){
			if (EStep1[1][j].indice != NULL){
				freelist(EStep1[1][j]);
			}
		}
		free(EStep1[1]);
		//used list deletion
				for (j=0;j<(1UL << l3);j++){
					if (EStep2[j].indice != NULL){
						freelist(EStep2[j]);
					}
				}
				free(EStep2);
	}
	else {
		//TODO l > 64 not implemented!
		printf("warning! l > 64 not implemented");
	}
	free(temp);
	free(temp2);
	free(indice);
	free(sums);
}

void sub_isd_report(unsigned long long cycles_per_iter) {
	(void) cycles_per_iter;
}

void sub_isd_free() {

}
