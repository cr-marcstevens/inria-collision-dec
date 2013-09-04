/**
 * \file BJMMtools.h
 * \brief Misc functions used by BJMM algorithm
 *
 * \author Mathieu Aria
 */

/**
 * \brief Linked list of sparse words + the coresponding sum
 */
typedef struct draw {
	word* sum;
	unsigned short* indice;
	struct draw* next;
} S ;

/**
 * \brief Compute next combination
 *
 * Can be called with 4 method (from 0 to 3) which share the matrix L with different patterns:
 * 
 *  - method 0 : [x------x||0------0] first half
 *  - method 1 : [x0x0x0x0x0x0x0x0x0] odd number
 *  - method 2 : [x-x||0-0||x-x||0-0] first and third quarter
 *  - method 3 : [x-x||0------0||x-x] first and final quarter
 *
 */
int next(word* L,unsigned short* indice,word* sums,unsigned int csize,unsigned int L_len, unsigned int w,int method);

int next2(word* L,unsigned short* indice,word* sums,unsigned int csize,unsigned int L_len, unsigned int w,int method);

void h1store(S* table,int shift,unsigned short* indice,word* sums,unsigned int w,unsigned int csize);

void fusionstore1(S* EStep1,word target,S* table,int shift1,int shift2,int shift3,unsigned short* indice,word* sums,unsigned int w,unsigned int csize);

void fusiongive1(S* answer,word target,S* table,int shift1,unsigned short* indice,word* sums,unsigned int w,unsigned int csize);

void FusionFilterStore64(S* AnswerList, S* StockedE,S* OnTheFlyE,word target,int shift1,int shift2,int eff_word_len,unsigned int l,unsigned int l2,unsigned int l3,unsigned int w,unsigned int w2,unsigned int csize);

void FusionFilterGive64(S** AnswerList, S* StockedE,S* OnTheFlyE,word target,int shift1,int shift2,unsigned int w,unsigned int w2,unsigned int csize);
#ifdef STAT
void FinalFusionFilter64(sw_list** AnswerList, S* StockedE,S* OnTheFlyE,word* Synd,int eff_word_len,unsigned int threshold,unsigned int l,unsigned int l2,unsigned int l3,unsigned int w,unsigned int w2,unsigned int csize, sw_list** observer);
#else
void FinalFusionFilter64(sw_list** AnswerList, S* StockedE,S* OnTheFlyE,word* Synd,int eff_word_len,unsigned int threshold,unsigned int l,unsigned int l2,unsigned int l3,unsigned int w,unsigned int w2,unsigned int csize);
#endif
void Sort(unsigned short* dest,unsigned short* s1,unsigned short* s2,unsigned short size1,unsigned short size2);

int SortFilter(unsigned short* dest,unsigned short* s1,unsigned short* s2,unsigned short size1,unsigned short size2,unsigned short targetsize);

void PrintDoubleStat(sw_list** AnswerList,unsigned short p);

void freelist(S* draw);
