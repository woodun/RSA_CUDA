
/******************************************************************************
 *
 *            (C) Copyright 2010 The Board of Trustees of the
 *                        University of Illinois
 *                         All Rights Reserved
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.cu"
#include <time.h>
#include "mpz.h"


int main (int argc, char *argv[])
{
	///////host memory
	long long int *clockTable_h;
	mpz_t *h_n;
	mpz_t *h_n_;
	mpz_t *h_r;
	mpz_t *h_r2;
	clockTable_h = (long long int*) malloc(10000 * sizeof(long long int));
	h_n = (mpz_t *) malloc (sizeof(mpz_t));
	h_n_ = (mpz_t*) malloc (sizeof(mpz_t));
	h_r = (mpz_t*) malloc (sizeof(mpz_t));
	mpz_init(h_n);
	mpz_init(h_n_);
	mpz_init(h_r);
	mpz_init(h_r2);

	///////get n
	char n_input[] = "38f6e8cfba55dd0e47";
	mpz_set_str_host(h_n, n_input, 1024)

	
	///////get n_
	char n__input[] = "2e8457440e0d93c489";
	mpz_set_str_host(h_n_, n__input, 1024)
	

	///////get r
	char r_input[] = "400000000000000000";
	mpz_set_str_host(h_r, r_input, 1024)


	///////get r2
	char r2_input[] = "3709d17d8f8686609f";
	mpz_set_str_host(h_r2, r2_input, 1024)


	///////get Messages
	int inputControl = 32;
	int mesSize = sizeof(mpz_t) * inputControl;
	mpz_t *myMes_h;
	myMes_h = (mpz_t*) malloc (mesSize);

	///////get Message1
	mpz_t *mes1;
	mes1 = (mpz_t*) malloc (sizeof(mpz_t));
	mpz_init(mes1);

	char mes1_input[] = ""; //input from pair storage
	mpz_set_str_host(mes1, mes1_input, 1024)
	for (int i=0; i<1; i++){
		mpz_init(&myMes_h[i]);
		mpz_set(&myMes_h[i], &mes1);
	}

	///////get Message2
	mpz_t *mes2;
	mes2 = (mpz_t*) malloc (sizeof(mpz_t));
	mpz_init(mes2);

	char mes2_input[] = ""; //input from pair storage
	mpz_set_str_host(mes2, mes2_input, 1024)
	for (int i=1; i<inputControl; i++){
		mpz_init(&myMes_h[i]);
		mpz_set(&myMes_h[i], &mes2);
	}

	mpz_t *myMes_d;
	cudaMalloc((mpz_t **) &myMes_d, mesSize);
	cudaMemcpy(myMes_d, myMes_h, mesSize, cudaMemcpyHostToDevice);


	///////get e
	char e_input[] = "1011";
	int e_bitsLength = (int)strlen(e_input);
	int* eBits = (int *) malloc(sizeof(int) * e_bitsLength);

	int* eBits_d;
	cudaMalloc((void **) &eBits_d, sizeof(int) * e_bitsLength);

	int e_iterator = e_bitsLength - 1;
	while ( e_iterator > 0){
        if( e_input == '1'){
            eBits[e_iterator] = 1;
        }
        else{
            eBits[e_iterator] = 0;
        }
        e_iterator--;
	}
	eBits[e_iterator] = 1;
	cudaMemcpy(eBits_d, eBits, sizeof(int) * e_bitsLength, cudaMemcpyHostToDevice);


	///////get d
	char d_input[] = "110111001111010000000011001001011000000011101110101111111100110011";
	int d_bitsLength = (int)strlen(d_input);
	int* dBits = (int *) malloc(sizeof(int) * d_bitsLength);

	int* eBits_d;
	cudaMalloc((void **) &eBits_d, sizeof(int) * d_bitsLength);

	int e_iterator = d_bitsLength - 1;
	while ( e_iterator > 0){
        if( e_input == '1'){
            dBits[e_iterator] = 1;
        }
        else{
            dBits[e_iterator] = 0;
        }
        e_iterator--;
	}
	dBits[e_iterator] = 1;
	cudaMemcpy(eBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);


	///////device memory 1
	mpz_t *tmp;
	mpz_t *tmp2;
	mpz_t *d_t;
	mpz_t *d_m;
	mpz_t *_a_mpz;
	mpz_t *_b_mpz;
	cudaMalloc((void **) &tmp, mesSize);
	cudaMalloc((void **) &tmp2, mesSize);
	cudaMalloc((void **) &d_t, mesSize);
	cudaMalloc((void **) &d_m, mesSize);
	cudaMalloc((void **) &_b_mpz, mesSize);
	cudaMalloc((void **) &_a_mpz, mesSize);


	///////device memory 2
	long long int *clockTable_d;
	mpz_t *d_n;
	mpz_t *d_n_;
	mpz_t *d_r;
	cudaMalloc((void **) &clockTable_d, 10000 * sizeof(long long int));
	cudaMalloc((void **) &d_n, sizeof(mpz_t));
	cudaMalloc((void **) &d_n_, sizeof(mpz_t));
	cudaMalloc((void **) &d_r, sizeof(mpz_t));


	cudaMemcpy(d_n, h_n, sizeof(mpz_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_n_, h_n_, sizeof(mpz_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, h_r, sizeof(mpz_t), cudaMemcpyHostToDevice);


	//MontSQMLadder<<<1, 32, 0>>>(myMes_d, _a_mpz, _b_mpz, tmp, tmp2, d_r, d_n, d_n_, eBits_d, e_bitsLength, clockTable_d, d_t, d_m);/////////////////////////////////////////kernel

	MontSQMLadder<<<1, 32, 0>>>(myMes_d, _a_mpz, _b_mpz, tmp, tmp2, d_r, d_n, d_n_, dBits_d, d_bitsLength, clockTable_d, d_t, d_m);/////////////////////////////////////////kernel
	cudaDeviceSynchronize();


	cudaMemcpy(clockTable_h, clockTable_d, 10000*sizeof(long long int), cudaMemcpyDeviceToHost);


	FILE *fp= fopen("oneCurve_allDeviant.txt", "w");
	for (int q = 0; q<10000; q++){
		fprintf(fp, "%lld\n", clockTable_h[q]);
	}
	fclose(fp);


	cudaFree(d_n);
	cudaFree(d_n_);
	cudaFree(d_r);
	cudaFree();//others
	cudaFree();
	cudaFree();
	cudaFree(myMes_d);
	free(myMes_h);
	free(clockTable_h);


    return 0;
}

