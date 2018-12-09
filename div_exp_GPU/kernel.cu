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

__device__ __host__ inline int CUDA_REDC(cuda_mpz_t* N, cuda_mpz_t* N_, cuda_mpz_t* T, cuda_mpz_t* tmp, cuda_mpz_t* t){//cuda_mpz_t* RMOD, int L, cuda_mpz_t* N, cuda_mpz_t* N_ should not be changed.

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
		return 1;
    }
	else{
		cuda_mpz_sub(tmp, t, N);
		cuda_mpz_set(tmp, t);
	    return 0;
	}
}

//__device__ __host__ inline cuda_mpz_t* CUDA_REDC(cuda_mpz_t* N, cuda_mpz_t* N_, cuda_mpz_t* T, cuda_mpz_t* tmp, cuda_mpz_t* t){//cuda_mpz_t* RMOD, int L, cuda_mpz_t* N, cuda_mpz_t* N_ should not be changed.
//
//	//m = ((T & R) * N_) & R
//	cuda_mpz_bitwise_truncate(t, T);
//	cuda_mpz_mult(tmp, N_, t);
//	cuda_mpz_bitwise_truncate_eq(tmp);
//
//	//t = (T + m*N) >> L
//	cuda_mpz_mult(t, tmp , N);
//	cuda_mpz_add(tmp, T, t);
//	cuda_mpz_bitwise_rshift(t, tmp);
//
//	unsigned carry = cuda_mpz_sub(tmp, t, N);
//	if (carry){
//		cuda_mpz_set(t, tmp);
//		return t;
//    }
//	else{
//		cuda_mpz_set(tmp, t);
//	    return t;
//	}
//}

__global__ void MontSQMLadder(cuda_mpz_t * mes1, cuda_mpz_t r2, cuda_mpz_t vn, cuda_mpz_t vn_, int* eBits, int eLength, long long int* divTable) {

	__shared__ cuda_mpz_t tmp[2]; //capacity will become a problem for shared memory with large keys
	__shared__ cuda_mpz_t tmp2[2];
	__shared__ cuda_mpz_t t[2];
	__shared__ cuda_mpz_t _x1[2];
	__shared__ cuda_mpz_t _x2[2];

	__shared__ int con1[2][1030];
	__shared__ int con2[2][1030];

	int j = threadIdx.x;
	int h = blockIdx.x * blockDim.x + threadIdx.x;

	if(j == 0){
		uint ret;
		asm("mov.u32 %0, %smid;" : "=r"(ret) );///make sure it's per SM
		printf("thread id: %d SM ID: %d\n", h, ret);//////////volta has 80 cores
	}

	///////////////////initialize
	cuda_mpz_init(&tmp[j]);
	cuda_mpz_init(&tmp2[j]);
	cuda_mpz_init(&t[j]);
	cuda_mpz_init(&_x1[j]);
	cuda_mpz_init(&_x2[j]);

	cuda_mpz_t* n = &vn;
	cuda_mpz_t* n_ = &vn_;

	cuda_mpz_set(&_x1[j], &mes1[h]);//next _x1 access will cause L1 miss if the L1 policy is write evict, same as using mutiple kernels.

	//_x1 = REDC(rmod,n,n_,mes*r2,l)
	cuda_mpz_mult(&tmp2[j], &_x1[j], &r2);
	con1[j][1024] = CUDA_REDC( n, n_, &tmp2[j], &tmp[j], &t[j]);
	cuda_mpz_set( &_x1[j], &t[j]);

	//x2 = _x1 * _x1
	cuda_mpz_mult(&tmp2[j], &_x1[j], &t[j]);
	//_x2 = CUDA_REDC(rmod,n,n_,_x2,l)
	con2[j][1024] = CUDA_REDC( n, n_, &tmp2[j], &tmp[j], &t[j]);
	cuda_mpz_set( &_x2[j], &t[j]);

//		if(j == 0){
//			printf("mes1: ");
//			cuda_mpz_print_str_device(&_x1[j]);
//			printf(" ");
//			cuda_mpz_print_str_device(&_x2[j]);
//			printf("\n");
//		}
//	for(int i = 1; i < eLength; ++i){///1 - 1023
//		if(eBits[i] == 0){
//			//x2 = _x1 * _x2
//			cuda_mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
//			//_x2 = CUDA_REDC(rmod,n,n_,_x2,l)
//			con1[j][i] = CUDA_REDC( n, n_, &tmp2[j], &tmp[j], &t[j]);
//			cuda_mpz_set( &_x2[j], &t[j]);
//			//_x1 = _x1 * _x1
//			cuda_mpz_set( &tmp[j], &_x1[j]);
//			cuda_mpz_mult(&tmp2[j], &_x1[j], &tmp[j]);
//			//_x1 = CUDA_REDC(rmod,n,n_,_x1,l)
//			con2[j][i] = CUDA_REDC( n, n_, &tmp2[j], &tmp[j], &t[j]);
//			cuda_mpz_set( &_x1[j], &t[j]);
//		} else {
//			//_x1 = _x1 * _x2
//			cuda_mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
//			//_x1 = CUDA_REDC(rmod,n,n_,_x1,l) #changes: more efficient
//			con1[j][i] = CUDA_REDC( n, n_, &tmp2[j], &tmp[j], &t[j]);
//			cuda_mpz_set( &_x1[j], &t[j]);
//			//_x2 = _x2 * _x2
//			cuda_mpz_set( &tmp[j], &_x2[j]);
//			cuda_mpz_mult(&tmp2[j], &_x2[j], &tmp[j]);
//			//_x2 = CUDA_REDC(rmod,n,n_,_x2,l) #changes: more efficient
//			con2[j][i] = CUDA_REDC( n, n_, &tmp2[j], &tmp[j], &t[j]);
//			cuda_mpz_set( &_x2[j], &t[j]);
//		}
//
	//_x1 = CUDA_REDC(rmod,n,n_,_x1,l)
	con1[j][1025] = CUDA_REDC( n, n_, &_x1[j], &tmp[j], &t[j]);
	cuda_mpz_set( &_x1[j], &t[j]);

	if( j == 0){
		int div_count = 0;
		for(int m = 1; m < 1025; m++){
			if(con1[0][m] != con1[1][m]){
				div_count++;
			}
			if(con2[0][m] != con2[1][m]){
				div_count++;
			}
		}
		if(con1[0][1025] != con1[1][1025]){
			div_count++;
		}
		divTable[h/2] = div_count;
	}
}







