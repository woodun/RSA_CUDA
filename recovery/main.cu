#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.cu"
#include <time.h>
#include "mpz.h"


long long unsigned time_diff(timespec start, timespec end){
	struct timespec temp;
	if ((end.tv_nsec - start.tv_nsec) < 0){
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}

	long long unsigned time_interval_ns = temp.tv_nsec;
	long long unsigned time_interval_s = temp.tv_sec;
	time_interval_s = time_interval_s * 1000000000;

	return time_interval_s + time_interval_ns;
}

//L1 enabled. (-keep -Xptxas -dlcm=ca --optimize 0)
int main (int argc, char *argv[]) {

	struct timespec ts1;
	clock_gettime(CLOCK_REALTIME, &ts1);

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

	mpz_t h_n;
	mpz_t h_n_;
	mpz_t h_r2;
	int rl = 70;

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
	long long unsigned mesSize = sizeof(mpz_t) * data_num;
//	printf("%llu %llu", data_num, mesSize);
	mpz_t *myMes1_h;
	myMes1_h = (mpz_t*) malloc (mesSize);

	for(long long unsigned i = 0; i < data_num; i++){
		mpz_init(&myMes1_h[i]);
	}

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
		mpz_set_str_host(&myMes1_h[line_num], line);
		line_num++;
		if(line_num == data_num){
			break;
		}
	}
	fclose(fp2);

	if (line)
	    free(line);

	mpz_t *myMes1_d;
	cudaMalloc((mpz_t **) &myMes1_d, mesSize);
	cudaMemcpy(myMes1_d, myMes1_h, mesSize, cudaMemcpyHostToDevice);

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
	char d_input[] = "1000100010110110111110111000110000000001011000001000011010101101000101";
	//char d_input[] = "101";

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
	unsigned varSize = sizeof(mpz_t) * thread_num;

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

	//int *divTable_d;///////div count
	//cudaMalloc((void **) &divTable_d, pairs * sizeof(int));///////div count
	//int *divTable_h;///////div count
	//divTable_h = (int*) malloc( pairs * sizeof(int));///////div count

	init<<<1, thread_num>>>(_x1_mpz, _x2_mpz, tmp, tmp2, d_t);
	cudaDeviceSynchronize();

	//MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, _x1_mpz, _x2_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t, divTable_d);///////div count kernel

	MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, _x1_mpz, _x2_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
	cudaDeviceSynchronize();

	cudaMemcpy(clockTable_h, clockTable_d, pairs * sizeof(long long int), cudaMemcpyDeviceToHost);
	//cudaMemcpy(divTable_h, divTable_d, pairs * sizeof(int), cudaMemcpyDeviceToHost);///////div count

	FILE *fp1= fopen(argv[3], "w");

	long long unsigned sum1 = 0;
	for (long long unsigned q = 256; q < pairs; q++){
		//fprintf(fp1, "%d ", divTable_h[q]);///////div count
		fprintf(fp1, "%lld\n", clockTable_h[q]);
		sum1 += clockTable_h[q];
	}
	sum1 = sum1 / (pairs - 256);

	printf("%lld", sum1);

	fclose(fp1);

	//cudaFree(divTable_d);///////div count
	//free(divTable_h);///////div count

	////////free device
	cudaFree(clockTable_d);
	cudaFree(eBits_d);
	cudaFree(dBits_d);
	cudaFree(myMes1_d);
	cudaFree(tmp);
	cudaFree(tmp2);
	cudaFree(d_t);
	cudaFree(_x1_mpz);
	cudaFree(_x2_mpz);

	////////free host
	free(clockTable_h);
	free(myMes1_h);
	free(eBits);
	free(dBits);

	/////////////////////////////////time
	struct timespec ts2;
	clock_gettime(CLOCK_REALTIME, &ts2);

	//printf("%llu ", time_diff(ts1, ts2));

    return 0;
}

