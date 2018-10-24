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


__device__ mpz_t* REDC(mpz_t* R,mpz_t* N,mpz_t* N_,mpz_t* T, mpz_t* tmp, mpz_t* tmp2, mpz_t* t, mpz_t* m){

	mpz_div(tmp, m, T, R);
	mpz_mult(tmp, m, N_);
	mpz_div(tmp2, m ,tmp, R);

	mpz_mult(t, m , N);
	mpz_add(tmp, T, t);
	mpz_div(t, tmp2,  tmp , R);

	if (mpz_gte(t , N)){
		mpz_sub(tmp, t, N);
		mpz_set(t, tmp);
		return t;

    }
	else{
		mpz_sub(tmp, t, N);
		mpz_set(tmp2, tmp);
	    return t;
	}
}

__global__ void MontSQMLadder(mpz_t * mes, mpz_t* _a, mpz_t* _b, mpz_t* tmp, mpz_t* tmp2, mpz_t* r, mpz_t* n, mpz_t* n_, int* eBits, int eLength, long long int* clockTable, mpz_t* t, mpz_t* m) {

	int j = blockIdx.x * blockDim.x + threadIdx.x;
	long long int t1, t2;
	mpz_init(&tmp[j]);
	mpz_init(&tmp2[j]);

	mpz_init(&_b[j]);
	mpz_init(&_a[j]);
	mpz_init(&t[j]);
	mpz_init(&m[j]);
	
	for(int a=0; a< 10000; a++){

		mpz_set_i(&_a[j], d_a);

		t1 = clock64(); //include as few lines as possible with the clocks

		mpz_mult(&tmp[j], &mes[j], r);

		mpz_div(&tmp2[j], &_b[j], &tmp[j], n);

		for(int i = 9; i>=0; i--){
			if(eBits[i] == 1){

				mpz_mult(&tmp[j], &_a[j], &_b[j]);
				mpz_set(&_a[j], &tmp[j]);

				mpz_set(&_a[j], REDC(r,n,n_,&_a[j],&tmp[j],&tmp2[j],&t[j],&m[j]));
			}

			mpz_mult(&tmp[j], &_b[j], &_b[j]);
			mpz_set(&_b[j], &tmp[j]);

			mpz_set(&_b[j], REDC(r,n,n_,&_b[j],&tmp[j],&tmp2[j],&t[j],&m[j]));
		}
		mpz_set(&_b[j], REDC(r,n,n_,&_a[j],&tmp[j], &tmp2[j],&t[j],&m[j]));

		t2 = clock64();
		
		if( j == 1){
			clockTable[a] = t2-t1;
		}
	}
}

