#include <stdio.h>
#include <stdlib.h>
#include "mpz.h"


__device__ mpz_cuda* REDC(int RL, mpz_cuda* N, mpz_cuda* N_, mpz_cuda* T, mpz_cuda* tmp, mpz_cuda* t){//mpz_cuda* RMOD, int L, mpz_cuda* N, mpz_cuda* N_ should not be changed.

	//m = ((T & R) * N_) & R
	mpz_bitwise_truncate(t, T, RL);
	mpz_mult(tmp, N_, t);
	mpz_bitwise_truncate_eq(tmp, RL);

	//t = (T + m*N) >> L
	mpz_mult(t, tmp , N);
	mpz_add(tmp, T, t);
	mpz_bitwise_rshift(t, tmp, RL);

	if (mpz_gte(t , N)){
		mpz_sub(tmp, t, N);
		mpz_set(t, tmp);
		return t;
    }
	else{
		mpz_sub(tmp, t, N);
		mpz_set(tmp, t);
	    return t;
	}
}

__global__ void MontSQMLadder(mpz_cuda * mes1, long long unsigned pairs, mpz_cuda* _x1, mpz_cuda* _x2, mpz_cuda* tmp, mpz_cuda* tmp2, int rl, mpz_cuda r2, mpz_cuda vn, mpz_cuda vn_, int* eBits, int eLength, long long int* clockTable, mpz_cuda* t) {

	__shared__ digit_t s_index[32];

	long long int t1, t2;

	long long int combo_num = 0;

	int k = blockIdx.x * blockDim.x + threadIdx.x;

	//to accelerate the experiment, we put all messages in one kernel launch. In the real case, each message causes one kernel launch.
	for(long long unsigned iter1 = 0; iter1 < pairs; iter1++){

		mpz_set(&_x1[k], &mes1[2 * iter1 + k]);//next _x1 access will cause L1 miss if the L1 policy is write evict, same as using mutiple kernels.
		s_index[k] = mpz_get_last_digit(&_x1[k]);//make a dependency to make sure previous store is finished.

		t1 = clock64();//beginning of necessary instructions within the kernel

		mpz_cuda* n = &vn;
		mpz_cuda* n_ = &vn_;
		int j = blockIdx.x * blockDim.x + threadIdx.x;

		//_x1 = REDC(rmod,n,n_,mes*r2,l)
		mpz_mult(&tmp2[j], &_x1[j], &r2);
		mpz_set( &_x1[j], REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]) );

		//x2 = _x1 * _x1
		mpz_mult(&tmp2[j], &_x1[j], &t[j]);
		//_x2 = REDC(rmod,n,n_,_x2,l)
		mpz_set( &_x2[j], REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]) );

//			if(j == 0){
//				mpz_print_str_device(&_x1[j]);
//				printf(" ");
//				mpz_print_str_device(&_x2[j]);
//				printf("\n");
//			}

		for(int i = eLength - 2; i >= 0; i--){

			if(eBits[i] == 0){
				//x2 = _x1 * _x2
				mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
				//_x2 = REDC(rmod,n,n_,_x2,l)
				mpz_set( &_x2[j], REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]) );
				//_x1 = _x1 * _x1
				mpz_set( &tmp[j], &_x1[j]);
				mpz_mult(&tmp2[j], &_x1[j], &tmp[j]);
				//_x1 = REDC(rmod,n,n_,_x1,l)
				mpz_set( &_x1[j], REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]) );
			} else {
				//_x1 = _x1 * _x2
				mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
				//_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
				mpz_set( &_x1[j], REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]) );
				//_x2 = _x2 * _x2
				mpz_set( &tmp[j], &_x2[j]);
				mpz_mult(&tmp2[j], &_x2[j], &tmp[j]);
				//_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
				mpz_set( &_x2[j], REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]) );
			}
		}

		//_x1 = REDC(rmod,n,n_,_x1,l)
		mpz_set( &_x1[j], REDC(rl, n, n_, &_x1[j], &tmp[j], &t[j]) );

		s_index[k] = mpz_get_last_digit(&_x1[k]);//make a dependency to make sure previous store is finished.

		t2 = clock64();//end of necessary kernel instructions

//		printf("combo_num: %lld, iter1: %u, iter2: %u\n", combo_num, iter1, iter2);

		if( j == 1){
			clockTable[iter1] = t2 - t1;
		}
	}
}

__global__ void init(mpz_cuda* _x1, mpz_cuda* _x2, mpz_cuda* tmp, mpz_cuda* tmp2, mpz_cuda* t){
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	mpz_init(&tmp[j]);////initial value not used
	mpz_init(&tmp2[j]);////initial value not used
	mpz_init(&_x1[j]);////initial value required
	mpz_init(&_x2[j]);////initial value required
	mpz_init(&t[j]);////initial value not used
}











