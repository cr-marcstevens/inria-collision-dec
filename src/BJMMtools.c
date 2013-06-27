/*
 * BJMMtools.c
 *
 *  Created on: 27 juin 2013
 *      Author: Mathieu ARIA
 */

#include "m4ri/m4ri.h"

/*
 * Cette fonction permet la construction d'un ensemble E de facon récursive à partir d'une pool de mot.
 * La methode indique le découpage de l'ensemble
 * method 1: |x-----x|0-----0|
 * method 2: |x0x0x0-------x0|
 * method 3: |x-x|0-0|x-x|0-0|
 * method 4: |x-x|0-----0|x-x|
 */
/*	void Ebuilder(word* pool,unsigned int poolsize,word* E,word current,unsigned int p,int method){
			int i;
			word* cible = E;
			switch(method){
			case 1:
				if (p==1){
					for(i=0;i<(poolsize-p+1);i++){
						(*(E+i)) = current ^ (*(pool+i));
						}
					}
				else{
					for(i=0;i<(poolsize-p+1);i++){
						Ebuilder(pool+i+1,poolsize-i-1,cible,(current ^ (*(pool+i))),p-1,method);
						cible += nCr(poolsize-i,p-1);
						}
					}
				}

			break;
			case 2: ;
			break;
			case 3: ;
			break;
			case 4: ;
			break;
			default: printf("building method unknown");
			}

	}
*/
