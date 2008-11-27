/*
 * Author: shadi lahham
 *
 * Main body that initiates chirality operations.
 *
 * deals with input output, and the main logic of the high level calculation
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> //for strcmp,strlen
#ifndef DDD
#include "permuter.h"
#endif
void nextPerm(int *perm, int size, int init, int num_used);
int temp;
int main() {
	int size = 13;
#ifndef RECURSIVE
	int i = 0;
	int c=0;
	permuter *p = createPermuter(size, 2);
	while (nextPermutation(p)) {		
//		for (i = 0; i < p->_size; i++) {
//			printf("%d ", p->_index[i] + 1);
//		}
		temp++;
//		printf("\n");
	}	
#else
	int* perm = (int *) malloc(sizeof(int) * size);
	int i = 0;
	for (i = 0; i < size; i++) {
		perm[i] = i;
	}
	nextPerm(perm, size, 0,0);	
#endif
	printf("done: %d\n", temp);
	return 0;
}

void nextPerm(int *perm, int size, int init,int num_used) {
	int i;
	if (num_used == size) {
		temp++;
//		for (i = 0; i < size; i++) {
//			printf ("%d ", perm[i] + 1);
//		}
//		printf("\n");
//		return;
	} 
	while(perm[init] != init && init < size) {
		init ++;
	}
	for (i = init; i < size; i++) {			
		if (perm[i] == i) {
			// space is clear
			int temp = perm[init];
			perm[init] = perm[i];
			perm[i] = temp;	
			nextPerm(perm, size, init + 1, num_used + ((i == init) ? 1 : 2));
			temp = perm[init];
			perm[init] = perm[i];
			perm[i] = temp;
		} else {
			nextPerm(perm, size, init + 1, num_used);
		}
	}
}
