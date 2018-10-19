/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpz.h"


__device__ mpz_t* REDC(mpz_t* R,mpz_t* N,mpz_t* N_,mpz_t* T, mpz_t* tmp, mpz_t* tmp2,mpz_t* t, mpz_t* m);


__global__ void reduction(mpz_t * mes, mpz_t* encryptMes, mpz_t* _a, mpz_t* _b, mpz_t* tmp, mpz_t* tmp2, mpz_t* r, mpz_t* n, mpz_t* n_, int* eBits, long long int* clockTable, mpz_t* t, mpz_t* m, int d_a) {


	int j = blockIdx.x * blockDim.x + threadIdx.x;
	long long int t1, t2;
	mpz_init(&tmp[j]);
	mpz_init(&tmp2[j]);

	mpz_init(&_b[j]);
	mpz_init(&_a[j]);
	mpz_init(&t[j]);
	mpz_init(&m[j]);

	
	if ( j < 32){
		
		for(int a=0; a< 10000; a++){
			t1 = clock64();
		
			mpz_set_i(&_a[j], d_a);			

			mpz_mult(&tmp[j], &mes[j], r);

			mpz_div(&tmp2[j], &_b[j], &tmp[j], n); // optimize
			//int _b = (mes[j]*r) % n;
		
		
			for(int i = 9; i>=0; i--){  // ebit traversal

				if(eBits[i] == 1){	// check whether ebit is one

					mpz_mult(&tmp[j], &_a[j], &_b[j]);
					mpz_set(&_a[j], &tmp[j]);
					//_a = _a * _b;
				
					mpz_set(&_a[j], REDC(r,n,n_,&_a[j],&tmp[j],&tmp2[j],&t[j],&m[j]));
					//_a =  REDC(r,n,n_,_a);
				}

				mpz_mult(&tmp[j], &_b[j], &_b[j]);
				mpz_set(&_b[j], &tmp[j]);
				//_b = _b * _b;
			
				mpz_set(&_b[j], REDC(r,n,n_,&_b[j],&tmp[j],&tmp2[j],&t[j],&m[j]));
				//_b = REDC(r,n,n_,_b);
			}
			mpz_set(&_b[j], REDC(r,n,n_,&_a[j],&tmp[j], &tmp2[j],&t[j],&m[j])); // final montgomery reduction, pass encrypted mes
		
			mpz_init(&encryptMes[j]); //redundant
			mpz_set(&encryptMes[j], &_b[j]); //redundant

			t2 = clock64();
			
			if( j == 1){
				clockTable[a] = t2-t1; // retrieve clock
			}
		}
	}
}

void Encrypt( mpz_t* mes, mpz_t* encrpytMes, int e, mpz_t* _a, mpz_t* _b, mpz_t* tmp, mpz_t* tmp2, mpz_t* r, mpz_t* n, mpz_t* n_, int* eBits, long long int* clockTable, mpz_t* t, mpz_t* m, int d_a){
//Encrypt (myMes_d, myMesEncrypted_d, e, _a_mpz, _b_mpz, tmp, tmp2, d_r, d_n, d_n_, eBits_d, clockTable_d, d_t, d_m, _a);
//Encrypt (myMes_d, myMesEncrypted_d, _a_mpz, _b_mpz, tmp, tmp2, d_r, d_n, d_n_, eBits_d, clockTable_d, d_t, d_m, _a);
	
	reduction<<<1, 32, 0>>>( mes, encrpytMes, _a, _b, tmp, tmp2, r, n, n_, eBits, clockTable, t, m, d_a);
	cudaDeviceSynchronize();
	
}


__device__ mpz_t* REDC(mpz_t* R,mpz_t* N,mpz_t* N_,mpz_t* T, mpz_t* tmp, mpz_t* tmp2, mpz_t* t, mpz_t* m){

	mpz_div(tmp, m, T, R);
	mpz_mult(tmp, m, N_); // k
	mpz_div(tmp2, m ,tmp, R);
	//int m = ((T % R) * N_) % R;
	
	mpz_mult(t, m , N);
	mpz_add(tmp, T, t);
	mpz_div(t, tmp2,  tmp , R);
	//int t = (T + m*N) / R;
	
	
	
	if (mpz_gte(t , N)){
		mpz_sub(tmp, t, N);
		mpz_set(t, tmp); //redundant
			
		return t;
		
    }
	else{
		mpz_sub(tmp, t, N); //redundant
		mpz_set(tmp2, tmp); //redundant
	    return t;
	}
}



__________
