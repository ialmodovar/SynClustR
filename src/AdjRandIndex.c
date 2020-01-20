/*************************************************
  @file: AdjRandIndex.c

Compute Adjusted Rand Index 
 
Author: 
Ranjan Maitra
Department of Statistics
Iowa State University
maitra@iastate.edu
 
Created in January 2018
**************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "array.h"
#define inf 1e+40;


/************************************************ 
This procedure computes Adjusted Rand index
  Parameters:
  N - number of observations
  TRUK - true number of clusters
  PREDK - estimated number of clusters
  trcl - true classification vector
  prcl - estimated classification vector
  Rand - Rand index
  adjRand - Adjusted Rand index
  Eindex - E index
**************************************************/

void RRand(int N, int TRUK, int PREDK, int *trcl, int *prcl, double *Rand, double *adjRand, double *F){
  
	// int i, j, n[TRUK][PREDK];
	int i, j;
	double sumtr[TRUK], sumpr[PREDK], sumprsq, sumtrsq, sumsq, discordant, sumtrprsq;
	double term1, term2, term3, W1, W2;
	double nij2sum, nidot2sum, ndotj2sum;

  /* Fixed by WCC. */
  int **n;
  n = (int **) malloc(TRUK * sizeof(int *));
  if(n == NULL) {
    perror("Memory allocation fails at n!\n");
  }
  for (i=0;i<TRUK;i++) {
    n[i] = (int *) malloc(PREDK * sizeof(int));
    if(n[i] == NULL) {
      perror("Memory allocation fails at n[i]!\n");
    }
    for (j=0;j<PREDK;j++) {
      n[i][j]=0;
    }
  }
  
	for (i = 0; i < N; i++){
		n[trcl[i]][prcl[i]] += 1;
	}

	sumtrsq = 0.;
	for (i = 0; i < TRUK; i++){
		sumtr[i] = 0.;
		for (j = 0; j < PREDK; j++){
			sumtr[i] += n[i][j];
		}
		sumtrsq += sumtr[i] * sumtr[i];
	}
  
	sumprsq = 0.;
	for (j = 0; j < PREDK; j++){
		sumpr[j] = 0.;
		for (i = 0; i < TRUK; i++){
			sumpr[j] += (double)n[i][j];
		}
		sumprsq += sumpr[j] * sumpr[j];
	}

	sumtrprsq = 0.;
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			sumtrprsq += sumtr[i] * sumtr[i] * sumpr[j] * sumpr[j];
		}
	}

//	(*Eindex) = sumtrprsq / (N * ((double)N - 1) + N * (double)N / (N - 1)) - (sumprsq + sumtrsq) / (N - 1);
//	(*Eindex) *= 2.;
//	(*Eindex) /= N * ((double)N - 1);
  
	sumsq = 0.;
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			sumsq += (double)n[i][j] * n[i][j];
		}
	}

	nij2sum = 0.;
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			nij2sum += (double)n[i][j] * (n[i][j] - 1) / 2.0;
		}
	}

	nidot2sum = 0.;
	for (i = 0; i < TRUK; i++){
		nidot2sum += (double)sumtr[i] * (sumtr[i] - 1) / 2.0;
	}

	ndotj2sum = 0.;
	for (i = 0; i < PREDK; i++){
		ndotj2sum += (double)sumpr[i] * (sumpr[i] - 1) / 2.0;
	}

	W1 = nij2sum / nidot2sum;
	W2 = nij2sum / ndotj2sum;
	
	(*F) = pow(W1 * W2, 0.5);
	
	discordant = 0.5 * (sumtrsq + sumprsq) - sumsq;

	(*Rand) = 1.0 - discordant / ((double)N * ((double)N - 1.) / 2.);

	term3 = nidot2sum * ndotj2sum / ((double)N * ((double)N - 1.) / 2.);

	term1 = nij2sum - term3;

	term2 = (nidot2sum + ndotj2sum) / 2 - term3;

	(*adjRand) = term1 / term2;

  /* Free 2D array pointers. */
  for (i=0;i<TRUK;i++){
    free(n[i]);
  }
  free(n);

}


void runAdjRand(int (*n), int (*K1), int (*K2), int *id1, int *id2, double (*Rand), double (*aRand), double (*F)){

	double Rand1, aRand1, F1;

	int n1, K11, K21;

	n1 = (*n);
	K11 = (*K1);
	K21 = (*K2);		

	Rand1 = (*Rand);
	aRand1 = (*aRand);
	F1 = (*F);
	
	RRand(n1, K11, K21, id1, id2, &Rand1, &aRand1, &F1);

	(*Rand) = Rand1;
	(*aRand) = aRand1;
	(*F) = F1;

}
