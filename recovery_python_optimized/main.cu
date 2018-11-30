#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.cu"
#include <time.h>
#include "cuda_mpz.h"


int main (int argc, char *argv[]) {

	///////input control
	if (argc < 2){
		exit(EXIT_FAILURE);
	}

	long x = strtol(argv[2], NULL, 10);
	long long unsigned pairs = x;
	unsigned thread_num = 2;
	long long unsigned data_num = pairs * 2;

	///////host memory
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc( pairs * sizeof(long long int));

	cuda_mpz_t h_n;
	cuda_mpz_t h_n_;
	cuda_mpz_t h_r2;

	///////get n
	char n_input[] = "00000038f6e8cfba55dd0e47";
	cuda_mpz_set_str_host(&h_n, n_input);
	
	///////get n_
	char n__input[] = "0000002e8457440e0d93c489";
	cuda_mpz_set_str_host(&h_n_, n__input);

	///////get r2
	char r2_input[] = "0000003709d17d8f8686609f";
	cuda_mpz_set_str_host(&h_r2, r2_input);

	///////get Messages
	long long unsigned mesSize = sizeof(cuda_mpz_t) * data_num;
	cuda_mpz_t *myMes1_h;
	myMes1_h = (cuda_mpz_t*) malloc (mesSize);

	///////get Message pairs
	char* line = NULL;
	size_t len = 0;

	FILE* fp2 = fopen(argv[1], "r");//input from pair storage
	if (fp2 == NULL){
	    exit(EXIT_FAILURE);
	}

	long long unsigned line_num = 0;
	while ((getline(&line, &len, fp2)) != -1) {
		line[strcspn(line, "\n")] = 0;
		cuda_mpz_set_str_host(&myMes1_h[line_num], line);
		line_num++;
		if(line_num == data_num){
			break;
		}
	}
	fclose(fp2);
	if (line)
	    free(line);

	cuda_mpz_t *myMes1_d;
	cudaMalloc((cuda_mpz_t **) &myMes1_d, mesSize);
	cudaMemcpy(myMes1_d, myMes1_h, mesSize, cudaMemcpyHostToDevice);

	///////get d
	char d_input[] = "1011011001001001010011110110010101010111001010110101111000111100001";

	int d_bitsLength = (int)strlen(d_input);
	int* dBits = (int *) malloc(sizeof(int) * d_bitsLength);
	int* dBits_d;
	cudaMalloc((void **) &dBits_d, sizeof(int) * d_bitsLength);

	int d_iterator = 0;
	while ( d_iterator < d_bitsLength){
		if( d_input[d_iterator] == '1'){//big endian
			dBits[d_iterator] = 1;
		}
		else{
			dBits[d_iterator] = 0;
		}
		d_iterator++;
	}
	cudaMemcpy(dBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);

	///////device memory
	long long int *clockTable_d;
	cudaMalloc((void **) &clockTable_d, 1 * sizeof(long long int));

	MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d);/////////////////////////////////////////kernel
	cudaDeviceSynchronize();

	cudaMemcpy(clockTable_h, clockTable_d, 1 * sizeof(long long int), cudaMemcpyDeviceToHost);
	long long int sum1 = clockTable_h[0];
	sum1 = sum1 / pairs;

	printf("%lld", sum1);

	////////free device
	cudaFree(clockTable_d);
	cudaFree(dBits_d);
	cudaFree(myMes1_d);

	////////free host
	free(clockTable_h);
	free(myMes1_h);
	free(dBits);

    return 0;
}

