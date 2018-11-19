#include <stdio.h>
#include <stdlib.h>
#include "mpz.h"


__device__ int REDC(int RL, mpz_t* N, mpz_t* N_, mpz_t* T, mpz_t* tmp, mpz_t* t){//mpz_t* RMOD, int L, mpz_t* N, mpz_t* N_ should not be changed.

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
		return 1;
    }
	else{
		mpz_sub(tmp, t, N);
		mpz_set(tmp, t);
	    return 0;
	}
}

__global__ void MontSQMLadder(mpz_t * mes1, long long unsigned pairs, mpz_t* _x1, mpz_t* _x2, mpz_t* tmp, mpz_t* tmp2, int rl, mpz_t r2, mpz_t vn, mpz_t vn_, int* eBits, int eLength, long long int* clockTable, mpz_t* t, int* divTable) {

	__shared__ digit_t s_index[32];
	__shared__ int con1[2][100];
	__shared__ int con2[2][100];

	long long int t1, t2;

	int k = blockIdx.x * blockDim.x + threadIdx.x;

	//to accelerate the experiment, we put all messages in one kernel launch. In the real case, each message causes one kernel launch.
	for(long long unsigned iter1 = 0; iter1 < pairs; iter1++){

		mpz_set(&_x1[k], &mes1[2 * iter1 + ((k + 30) / 31) ]);//next _x1 access will cause L1 miss if the L1 policy is write evict, same as using mutiple kernels.
		s_index[k] = mpz_get_last_digit(&_x1[k]);//make a dependency to make sure previous store is finished.

		t1 = clock64();//beginning of necessary instructions within the kernel

		mpz_t* n = &vn;
		mpz_t* n_ = &vn_;
		int j = blockIdx.x * blockDim.x + threadIdx.x;

		//_x1 = REDC(rmod,n,n_,mes*r2,l)
		mpz_mult(&tmp2[j], &_x1[j], &r2);
		con1[j][70] = REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]);
		mpz_set( &_x1[j], &t[j]);

		//x2 = _x1 * _x1
		mpz_mult(&tmp2[j], &_x1[j], &t[j]);
		//_x2 = REDC(rmod,n,n_,_x2,l)
		con2[j][70] = REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]);
		mpz_set( &_x2[j], &t[j]);

//			if(j == 0){
//				mpz_print_str_device(&_x1[j]);
//				printf(" ");
//				mpz_print_str_device(&_x2[j]);
//				printf("\n");
//			}

//		if(j == 0){
//			printf("(%)",con1[j][70], con1[j][70], con1[j][70], con1[j][70]);
//		}

		for(int i = eLength - 2; i >= 0; i--){ //0 - 65 = 67 - 1 //0 - 68 = 70 - 1

			if(eBits[i] == 0){
				//x2 = _x1 * _x2
				mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
				//_x2 = REDC(rmod,n,n_,_x2,l)
				con1[j][i] = REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]);
				mpz_set( &_x2[j], &t[j]);
				//_x1 = _x1 * _x1
				mpz_set( &tmp[j], &_x1[j]);
				mpz_mult(&tmp2[j], &_x1[j], &tmp[j]);
				//_x1 = REDC(rmod,n,n_,_x1,l)
				con2[j][i] = REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]);
				mpz_set( &_x1[j], &t[j]);
			} else {
				//_x1 = _x1 * _x2
				mpz_mult(&tmp2[j], &_x1[j], &_x2[j]);
				//_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
				con1[j][i] = REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]);
				mpz_set( &_x1[j], &t[j]);
				//_x2 = _x2 * _x2
				mpz_set( &tmp[j], &_x2[j]);
				mpz_mult(&tmp2[j], &_x2[j], &tmp[j]);
				//_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
				con2[j][i] = REDC(rl, n, n_, &tmp2[j], &tmp[j], &t[j]);
				mpz_set( &_x2[j], &t[j]);
			}
		}

		//_x1 = REDC(rmod,n,n_,_x1,l)
		con1[j][69] = REDC(rl, n, n_, &_x1[j], &tmp[j], &t[j]);
		mpz_set( &_x1[j], &t[j]);

		s_index[k] = mpz_get_last_digit(&_x1[k]);//make a dependency to make sure previous store is finished.

		t2 = clock64();//end of necessary kernel instructions

		if( j == 1){
			clockTable[iter1] = t2 - t1;

			con2[0][69] = 0;
			con2[1][69] = 0;
			int div_count = 0;
			for(int m = 0; m < 71; m++){
				if(con1[0][m] != con1[1][m]){
					div_count++;
				}
				if(con2[0][m] != con2[1][m]){
					div_count++;
				}
			}
			divTable[iter1] = div_count;////////////////////////////////////////////////////////////////todo: where is wrong with this number?
		}
	}
}

__global__ void init(mpz_t* _x1, mpz_t* _x2, mpz_t* tmp, mpz_t* tmp2, mpz_t* t){
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	mpz_init(&tmp[j]);////initial value not used
	mpz_init(&tmp2[j]);////initial value not used
	mpz_init(&_x1[j]);////initial value required
	mpz_init(&_x2[j]);////initial value required
	mpz_init(&t[j]);////initial value not used
}











