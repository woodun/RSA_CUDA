
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
	clockTable_h = (long long int*) malloc(10000 * sizeof(long long int));

	mpz_t h_n;
	mpz_t h_n_;
	mpz_t h_r2;
	int rl = 70;
	unsigned mes_size = 1;
	int inputControl = 32;

	mpz_init(&h_n);
	mpz_init(&h_n_);
	mpz_init(&h_r2);

	///////get n
	char n_input[] = "00000038f6e8cfba55dd0e47";
	mpz_set_str_host(&h_n, n_input);
	
	///////get n_
	char n__input[] = "0000000038f6e8cfba55dd0e47";
	mpz_set_str_host(&h_n_, n__input);

	///////get r2
	char r2_input[] = "0000000038f6e8cfba55dd0e47";
	mpz_set_str_host(&h_r2, r2_input);

	///////get Messages
	int mesSize = sizeof(mpz_t) * inputControl;
	mpz_t *myMes_h;
	myMes_h = (mpz_t*) malloc (mesSize);

	///////get Message1
	char mes1_input[] = "00012345"; //input from pair storage
	mpz_set_str_host(&myMes_h[0], mes1_input);

	///////get Message2
	char mes2_input[] = "00067890"; //input from pair storage
	mpz_set_str_host(&myMes_h[1], mes2_input);

	//for (int i=0; i<1; i++){
	//	mpz_init(&myMes_h[i]);
	//	mpz_set(&myMes_h[i], &mes1);
	//}

	//char test_str[1024];
	//printf("%s\n", mpz_get_str(&h_n, test_str, 1024));
	//exit(0);

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
        if( e_input[e_bitsLength - 1 - e_iterator] == '1'){
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

	int* dBits_d;
	cudaMalloc((void **) &dBits_d, sizeof(int) * d_bitsLength);

	int d_iterator = d_bitsLength - 1;
	while ( d_iterator > 0){

		printf("%c", dBits[d_input[e_bitsLength - 1 - d_iterator]]);

        if( d_input[e_bitsLength - 1 - d_iterator] == '1'){
            dBits[d_iterator] = 1;
            printf("1");
        }
        else{
            dBits[d_iterator] = 0;
            printf("0");
        }
        d_iterator--;
	}
	dBits[d_iterator] = 1;
	cudaMemcpy(dBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);

	printf("\n%d\n", d_bitsLength);
	for (int i = 0; i < d_bitsLength; i++){
		printf("%d", dBits[i]);
	}
	exit(0);

	///////device memory
	long long int *clockTable_d;
	mpz_t *tmp;
	mpz_t *tmp2;
	mpz_t *d_t;
	mpz_t *_x1_mpz;
	mpz_t *_x2_mpz;
	cudaMalloc((void **) &clockTable_d, 10000 * sizeof(long long int));
	cudaMalloc((void **) &tmp, mesSize);
	cudaMalloc((void **) &tmp2, mesSize);
	cudaMalloc((void **) &d_t, mesSize);
	cudaMalloc((void **) &_x2_mpz, mesSize);
	cudaMalloc((void **) &_x1_mpz, mesSize);

	//init(mpz_t* _x1, mpz_t* _x2, mpz_t* tmp, mpz_t* tmp2, mpz_t* t){
	init<<<1, inputControl>>>(_x1_mpz, _x2_mpz, tmp, tmp2, d_t);
	cudaDeviceSynchronize();

	//MontSQMLadder<<<1, inputControl>>>(myMes_d, mes_size, _x1_mpz, _x2_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, eBits_d, e_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
	//cudaDeviceSynchronize();

	//MontSQMLadder(mpz_t * mes, unsigned mes_size, mpz_t* _x1, mpz_t* _x2, mpz_t* tmp, mpz_t* tmp2, int rl, mpz_t r2, mpz_t vn, mpz_t vn_, int* eBits, int eLength, long long int* clockTable, mpz_t* t)
	MontSQMLadder<<<1, inputControl>>>(myMes_d, mes_size, _x1_mpz, _x2_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
	cudaDeviceSynchronize();

	cudaMemcpy(clockTable_h, clockTable_d, 10000*sizeof(long long int), cudaMemcpyDeviceToHost);

	FILE *fp= fopen("oneCurve_allDeviant.txt", "w");
	for (int q = 0; q<10000; q++){
		fprintf(fp, "%lld\n", clockTable_h[q]);
	}
	fclose(fp);

	////////free device
	cudaFree(clockTable_d);
	cudaFree(eBits_d);
	cudaFree(dBits_d);
	cudaFree(myMes_d);
	cudaFree(tmp);
	cudaFree(tmp2);
	cudaFree(d_t);
	cudaFree(_x1_mpz);
	cudaFree(_x2_mpz);

	////////free host
	free(clockTable_h);
	free(myMes_h);
	free(eBits);
	free(dBits);

    return 0;
}

