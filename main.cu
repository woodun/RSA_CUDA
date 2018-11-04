
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

	FILE* fp2 = fopen(argv[2], "r");
	if (fp2 == NULL)
	    exit(EXIT_FAILURE);

	char* line = NULL;
	size_t len = 0;
	mpz_t test;
	mpz_init(&test);
	while ((getline(&line, &len, fp2)) != -1) {
		printf("%s", line);
		mpz_set_str_host(&test, line);
		char test_str[1024];
		printf("%s\n", mpz_get_str(&test, test_str, 1024));
	}
	fclose(fp2);
	if (line)
	    free(line);

	exit(0);//./main asd.txt 1branchcombo0000.txt

	///////host memory
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc(10000 * sizeof(long long int));

	mpz_t h_n;
	mpz_t h_n_;
	mpz_t h_r2;
	int rl = 70;
	unsigned pairs = 1002001;
	unsigned inputControl = 2;
	unsigned mes_number = pairs * inputControl;

	mpz_init(&h_n);
	mpz_init(&h_n_);
	mpz_init(&h_r2);

	///////get n
	char n_input[] = "00000038f6e8cfba55dd0e47";
	mpz_set_str_host(&h_n, n_input);
	
	///////get n_
	char n__input[] = "0000002e8457440e0d93c489";
	mpz_set_str_host(&h_n_, n__input);

	///////get r2
	char r2_input[] = "0000003709d17d8f8686609f";
	mpz_set_str_host(&h_r2, r2_input);

	///////get Messages
	unsigned mesSize = sizeof(mpz_t) * mes_number;
	mpz_t *myMes_h;
	myMes_h = (mpz_t*) malloc (mesSize);

	for(int i = 0; i < mes_number; i++){
		mpz_init(&myMes_h[i]);
	}

	///////get Message1
	char mes1_input[] = "00000000000123456789"; //input from pair storage
	mpz_set_str_host(&myMes_h[0], mes1_input);

	///////get Message2
	char mes2_input[] = "00000000000987654321"; //input from pair storage
	mpz_set_str_host(&myMes_h[1], mes2_input);

//	//debug
//	char test_str[1024];
//	printf("%s\n", mpz_get_str(&h_n, test_str, 1024));
//	printf("%s\n", mpz_get_str(&h_n_, test_str, 1024));
//	printf("%s\n", mpz_get_str(&h_r2, test_str, 1024));

	mpz_t *myMes_d;
	cudaMalloc((mpz_t **) &myMes_d, mesSize);
	cudaMemcpy(myMes_d, myMes_h, mesSize, cudaMemcpyHostToDevice);

	///////get e
	char e_input[] = "101";
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
	char d_input[] = "1011011001001001010011110110010101010111001010110101111000111100001";

	int d_bitsLength = (int)strlen(d_input);

	int* dBits = (int *) malloc(sizeof(int) * d_bitsLength);

	int* dBits_d;
	cudaMalloc((void **) &dBits_d, sizeof(int) * d_bitsLength);

	int d_iterator = d_bitsLength - 1;
	while ( d_iterator > 0){
        if( d_input[d_bitsLength - 1 - d_iterator] == '1'){
            dBits[d_iterator] = 1;
        }
        else{
            dBits[d_iterator] = 0;
        }
        d_iterator--;
	}
	dBits[d_iterator] = 1;
	cudaMemcpy(dBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);

	///////device memory
	unsigned varSize = sizeof(mpz_t) * inputControl;

	long long int *clockTable_d;
	mpz_t *tmp;
	mpz_t *tmp2;
	mpz_t *d_t;
	mpz_t *_x1_mpz;
	mpz_t *_x2_mpz;
	cudaMalloc((void **) &clockTable_d, pairs * sizeof(long long int));
	cudaMalloc((void **) &tmp, varSize);
	cudaMalloc((void **) &tmp2, varSize);
	cudaMalloc((void **) &d_t, varSize);
	cudaMalloc((void **) &_x2_mpz, varSize);
	cudaMalloc((void **) &_x1_mpz, varSize);

	init<<<1, inputControl>>>(_x1_mpz, _x2_mpz, tmp, tmp2, d_t);
	cudaDeviceSynchronize();

//	printf("x1: %s\n", mpz_get_str(&myMes_h[0], test_str, 1024));
//	printf("x2: %s\n", mpz_get_str(&myMes_h[1], test_str, 1024));
//
//	MontSQMLadder<<<1, inputControl>>>(myMes_d, pairs, _x1_mpz, _x2_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, eBits_d, e_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
//	cudaDeviceSynchronize();
//
//	cudaMemcpy(myMes_d, _x1_mpz, mesSize, cudaMemcpyDeviceToDevice);
//	cudaMemcpy(myMes_h, _x1_mpz, mesSize, cudaMemcpyDeviceToHost);
//
//	printf("x1: %s\n", mpz_get_str(&myMes_h[0], test_str, 1024));
//	printf("x2: %s\n", mpz_get_str(&myMes_h[1], test_str, 1024));

	MontSQMLadder<<<1, inputControl>>>(myMes_d, pairs, _x1_mpz, _x2_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
	cudaDeviceSynchronize();

//	cudaMemcpy(myMes_h, _x1_mpz, mesSize, cudaMemcpyDeviceToHost);
//
//	printf("x1: %s\n", mpz_get_str(&myMes_h[0], test_str, 1024));
//	printf("x2: %s\n", mpz_get_str(&myMes_h[1], test_str, 1024));

	cudaMemcpy(clockTable_h, clockTable_d, pairs * sizeof(long long int), cudaMemcpyDeviceToHost);

	FILE *fp= fopen(argv[1], "w");
	for (int q = 0; q < pairs; q++){
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

