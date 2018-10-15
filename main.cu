
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


struct egcdTuple {
	int g;
	int y;
	int x;
};///gcd a/gcd b/gcd

egcdTuple egcd(int a, int b, struct egcdTuple c);
int modinv(int a, int m, struct egcdTuple c);
mpz_t* ABC(mpz_t* R,mpz_t* N,mpz_t* N_,mpz_t* T, mpz_t* tmp, mpz_t* tmp2, mpz_t* d_t, mpz_t* d_m);



int main (int argc, char *argv[])
{

    // Initialize host variables ----------------------------------------------
    int inputControl = 32;
	mpz_t *myMes_h;
	mpz_t *myMes_d;
	mpz_t *myMesEncrypted_h;
	mpz_t *myMesEncrypted_d;
	mpz_t *_a_mpz;
	mpz_t *_b_mpz;
	
	
	int mesSize = sizeof(mpz_t) * inputControl;
	long long int *clockTable_h;
	long long int *clockTable_d;
	long long int *clockTableT1_h;
	long long int *clockTableT1_d;
	long long int *clockTableT2_h;
	long long int *clockTableT2_d;
	int clockSize = sizeof(long long int) * inputControl;
	
	int p = 107;
	int q = 83;
	int n = p*q; // = 8881
	int e = 3;
	int temp = n;
	int length = 1;
	
	while (temp > 1){
		temp = temp / 2;
		length += 1;

	}
	
	// a generation
	int r = pow(2,length);
	int n_;
	struct egcdTuple modInvVariable;
	n_ = - modinv(n,r, modInvVariable) % r;
	int _a = (1*r) % n;	
	
	// e generation
	int bitsLength = 10;
	int* eBits = (int *) malloc(sizeof(int) * bitsLength);
	int e_iterator = bitsLength - 1;
	int e_temp = e;
	int* eBits_d;
	cudaMalloc((void **) &eBits_d, sizeof(int) * bitsLength);

	
	for (int i=0; i<10; i++){
        eBits[i] = 0;
	}
	

	while ( e_temp > 1){
        if( e_temp%2 == 1){
            eBits[e_iterator] = 1;
        }
        else{
            eBits[e_iterator] = 0;
        }
        e_iterator--;
        e_temp = e_temp / 2;
	}
	eBits[e_iterator] = 1;
	
	cudaMemcpy(eBits_d, eBits, sizeof(int) * bitsLength, cudaMemcpyHostToDevice);	
	cudaDeviceSynchronize();
	

	// Mes transfer setup
	myMes_h = (mpz_t*) malloc (mesSize);
	myMesEncrypted_h = (mpz_t*) malloc (mesSize);
	clockTable_h = (long long int*) malloc(10000 * sizeof(long long int));
	clockTableT1_h = (long long int*) malloc(clockSize);
	clockTableT2_h = (long long int*) malloc(clockSize);
	
	
	
	mpz_init(myMes_h);//redundant?
	mpz_init(myMesEncrypted_h);//redundant?

	cudaMalloc((mpz_t **) &myMes_d, mesSize);
	cudaMalloc((mpz_t **) &myMesEncrypted_d, mesSize);
	cudaMalloc((void **) &clockTable_d, 10000 * sizeof(long long int));
	cudaMalloc((void **) &clockTableT1_d, clockSize);
	cudaMalloc((void **) &clockTableT2_d, clockSize);

	//cudaDeviceSynchronize();
	
	
	//initializeArray(myMes_d, myMesEncrypted_d);
	//cudaDeviceSynchronize();
	
	
	//&myMes_h[1] = (mpz_t*) malloc(sizeof(_MPZ));
	
	for (int i=0; i<inputControl; i++){
		mpz_init(&myMes_h[i]);
		mpz_init(&myMesEncrypted_h[i]);
		mpz_set_i(&myMes_h[i], 123 + i);

	}

	
	
	
	cudaMemcpy(myMes_d, myMes_h, mesSize, cudaMemcpyHostToDevice);	
	
	//std::fstream fs;
	//fs.open ("test_uniform.txt", std::fstream::in | std::fstream::out | std::ofstream::trunc| std::fstream::app );
	FILE *fp= fopen("oneCurve_allDeviant.txt", "w");
	
	
	
	
	mpz_t *tmp;
	mpz_t *tmp2;
	mpz_t *d_n;
	mpz_t *d_n_;
	mpz_t *h_n;
	mpz_t *h_n_;
	mpz_t *d_r;
	mpz_t *h_r;
	mpz_t *d_t;
	mpz_t *d_m;
	cudaMalloc((void **) &d_n, sizeof(mpz_t));
	cudaMalloc((void **) &d_n_, sizeof(mpz_t));
	cudaMalloc((void **) &d_r, sizeof(mpz_t));
	cudaMalloc((void **) &d_t, mesSize);
	cudaMalloc((void **) &d_m, mesSize);
	cudaMalloc((void **) &tmp, mesSize);
	cudaMalloc((void **) &tmp2, mesSize);
	cudaMalloc((void **) &_b_mpz, mesSize);
	cudaMalloc((void **) &_a_mpz, mesSize);

	h_n = (mpz_t *) malloc (sizeof(mpz_t));
	h_n_ = (mpz_t*) malloc (sizeof(mpz_t));
	h_r = (mpz_t*) malloc (sizeof(mpz_t));
	
	mpz_init(h_n);
	mpz_init(h_n_);
	mpz_init(h_r);
	
	mpz_set_i(h_n, n);
	mpz_set_i(h_n_, n_);
	mpz_set_i(h_r, r);
	
	cudaMemcpy(d_n, h_n, sizeof(mpz_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_n_, h_n_, sizeof(mpz_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, h_r, sizeof(mpz_t), cudaMemcpyHostToDevice);

	

	Encrypt (myMes_d, myMesEncrypted_d, e, _a_mpz, _b_mpz, tmp, tmp2, d_r, d_n, d_n_, eBits_d, clockTable_d, d_t, d_m, _a);
	

	cudaMemcpy(clockTable_h, clockTable_d, 10000*sizeof(long long int), cudaMemcpyDeviceToHost);

	//fs  <<  clockTable_h[0] << "\n";
	for (int q = 0; q<10000; q++){
		fprintf(fp, "%lld\n", clockTable_h[q]);
	}
	/*
	for (int z = 0; z<inputControl; z++){
		printf("\n%d th input is \n", z);
		mpz_print(&myMes_h[z]);
		mpz_print(&myMesEncrypted_h[z]);
		//printf("clock table is %lld \n", clockTable_h[z]);
		//printf("T1 is %lld and T2 is %lld", clockTableT1_h[z], clockTableT2_h[z]);
	}*/
	
	
	//fs.close();
	fclose(fp);

	cudaFree(myMesEncrypted_d);
	cudaFree(myMes_d);
	free(myMes_h);
	free(myMesEncrypted_h);
	free(clockTable_h);

	

    return 0;

}

egcdTuple egcd(int a, int b, struct egcdTuple c){
    if (a == 0){
		return {b,0,1};
    }
	else{

        c = egcd( b%a, a, c);
		return {c.g, c.x - (b/a) * c.y, c.y};
	}
}


int modinv(int a, int m, struct egcdTuple c){
    c = egcd(a, m, c);
    int temp_x, temp_y;
    temp_x = c.x;
    temp_y = c.y;
    c.x = temp_y;
    c.y = temp_x;

    if (c.g != 1) {
    	return 0;//raise Exception('modular inverse does not exist');
	}
    else{

        return (c.x % m);
	}
}
