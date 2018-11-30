#include <stdio.h>
#include <stdlib.h>
#include "cuda_mpz.h"


__device__ __host__ inline cuda_mpz_t* REDC(cuda_mpz_t* N, cuda_mpz_t* N_, cuda_mpz_t* T, cuda_mpz_t* tmp, cuda_mpz_t* t){//cuda_mpz_t* RMOD, int L, cuda_mpz_t* N, cuda_mpz_t* N_ should not be changed.

	//m = ((T & R) * N_) & R
	cuda_mpz_bitwise_truncate(t, T);
	cuda_mpz_mult(tmp, N_, t);
	cuda_mpz_bitwise_truncate_eq(tmp);

	//t = (T + m*N) >> L
	cuda_mpz_mult(t, tmp , N);
	cuda_mpz_add(tmp, T, t);
	cuda_mpz_bitwise_rshift(t, tmp);

	if (cuda_mpz_gte(t , N)){
		cuda_mpz_sub(tmp, t, N);
		cuda_mpz_set(t, tmp);
		return t;
    }
	else{
		cuda_mpz_sub(tmp, t, N);
		cuda_mpz_set(tmp, t);
	    return t;
	}
}

__global__ void MontSQMLadder(cuda_mpz_t * mes1, long long unsigned pairs, cuda_mpz_t r2, cuda_mpz_t vn, cuda_mpz_t vn_, int* eBits, int eLength, long long int* clockTable) {

	__shared__ cuda_mpz_t tmp[2];
	__shared__ cuda_mpz_t tmp2[2];
	__shared__ cuda_mpz_t t[2];
	__shared__ cuda_mpz_t _x1[2];
	__shared__ cuda_mpz_t _x2[2];

	__shared__ digit_t s_index[32];

	long long int t1, t2;
	long long int sum1 = 0;

	int k = blockIdx.x * blockDim.x + threadIdx.x;

	///////////////////initialize
	cuda_mpz_init(&tmp[k]);
	cuda_mpz_init(&tmp2[k]);
	cuda_mpz_init(&t[k]);
	cuda_mpz_init(&_x1[k]);
	cuda_mpz_init(&_x2[k]);

	cuda_mpz_t* n = &vn;
	cuda_mpz_t* n_ = &vn_;

	//to accelerate the experiment, we put all messages in one kernel launch. In the real case, each message causes one kernel launch.
	for(unsigned iter1 = 0; iter1 < pairs * 4; ++iter1){

		s_index[k] = iter1;//make a dependency to make sure previous store is finished.
		t1 = clock64();//beginning of necessary instructions within the kernel

		cuda_mpz_set(&_x1[k], &mes1[2 * iter1 + k]);//next _x1 access will cause L1 miss if the L1 policy is write evict, same as using mutiple kernels.

		int j = blockIdx.x * blockDim.x + threadIdx.x;

		//_x1 = REDC(rmod,n,n_,mes*r2,l)
		cuda_mpz_mult(&tmp2[j], &_x1[j], &r2);
		cuda_mpz_set( &_x1[j], REDC(n, n_, &tmp2[j], &tmp[j], &t[j]) );
		//x2 = _x1 * _x1
		cuda_mpz_mult(&tmp2[j], &_x1[j], &t[j]);
		//_x2 = REDC(rmod,n,n_,_x2,l)
		cuda_mpz_set( &_x2[j], REDC(n, n_, &tmp2[j], &tmp[j], &t[j]) );

//		if(j == 0){
//			cuda_mpz_print_str_device(&_x1[j]);
//			printf(" ");
//			cuda_mpz_print_str_device(&_x2[j]);
//			printf("\n");
//		}

		for(int i = 1; i < eLength; ++i){

			if(eBits[i] == 0){
				//x2 = _x1 * _x2
				cuda_mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
				//_x2 = REDC(rmod,n,n_,_x2,l)
				cuda_mpz_set( &_x2[j], REDC(n, n_, &tmp2[j], &tmp[j], &t[j]) );
				//_x1 = _x1 * _x1
				cuda_mpz_set( &tmp[j], &_x1[j]);
				cuda_mpz_mult(&tmp2[j], &_x1[j], &tmp[j]);
				//_x1 = REDC(rmod,n,n_,_x1,l)
				cuda_mpz_set( &_x1[j], REDC(n, n_, &tmp2[j], &tmp[j], &t[j]) );
			} else {
				//_x1 = _x1 * _x2
				cuda_mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
				//_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
				cuda_mpz_set( &_x1[j], REDC(n, n_, &tmp2[j], &tmp[j], &t[j]) );
				//_x2 = _x2 * _x2
				cuda_mpz_set( &tmp[j], &_x2[j]);
				cuda_mpz_mult(&tmp2[j], &_x2[j], &tmp[j]);
				//_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
				cuda_mpz_set( &_x2[j], REDC(n, n_, &tmp2[j], &tmp[j], &t[j]) );
			}
		}

		//_x1 = REDC(rmod,n,n_,_x1,l)
		cuda_mpz_set( &_x1[j], REDC(n, n_, &_x1[j], &tmp[j], &t[j]) );

		s_index[k] = cuda_mpz_get_last_digit(&_x1[k]);//make a dependency to make sure previous store is finished.
		t2 = clock64();//end of necessary kernel instructions
//		printf("combo_num: %lld, iter1: %u, iter2: %u\n", combo_num, iter1, iter2);

		sum1 +=  t2 - t1;
		if (iter1 == pairs - 1){
			if( k == 0){
				clockTable[0] = sum1;
			}
		}
		if (iter1 == 2 * pairs - 1){
			if( k == 0){
				clockTable[1] = sum1;
			}
		}
	}
}







