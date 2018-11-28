#include <stdio.h>
#include <stdlib.h>
//#include "mpz.h"
#include "cuda_mpz.h"

__host__ mpz_t* REDC(int RL, mpz_t* N, mpz_t* N_, mpz_t* T, mpz_t* tmp, mpz_t* t){//mpz_t* RMOD, int L, mpz_t* N, mpz_t* N_ should not be changed.

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

void MontSQMLadder(mpz_t * mes1, long long unsigned pairs, mpz_t* _x1, mpz_t* _x2, mpz_t* tmp, mpz_t* tmp2, int rl, mpz_t r2, mpz_t vn, mpz_t vn_, int* eBits, int eLength, mpz_t* t) {

	//to accelerate the experiment, we put all messages in one kernel launch. In the real case, each message causes one kernel launch.
	for(long long unsigned iter1 = 0; iter1 < pairs; iter1++){

		mpz_set(&_x1[0], &mes1[0]);//next _x1 access will cause L1 miss if the L1 policy is write evict, same as using mutiple kernels.

		mpz_t* n = &vn;
		mpz_t* n_ = &vn_;

//		if(j == 1){
//			printf("mes2: ");
//			mpz_print_str_device(&_x1[0]);
//			printf("\n");
//		}

		//_x1 = REDC(rmod,n,n_,mes*r2,l)
		mpz_mult(&tmp2[0], &_x1[0], &r2);
		mpz_set( &_x1[0], REDC(rl, n, n_, &tmp2[0], &tmp[0], &t[0]) );

		//x2 = _x1 * _x1
		mpz_mult(&tmp2[0], &_x1[0], &t[0]);
		//_x2 = REDC(rmod,n,n_,_x2,l)
		mpz_set( &_x2[0], REDC(rl, n, n_, &tmp2[0], &tmp[0], &t[0]) );

//		if(j == 0){
//			mpz_print_str_device(&_x1[0]);
//			printf(" ");
//			mpz_print_str_device(&_x2[0]);
//			printf("\n");
//		}

//		if(j == 0){
//			printf("mes1: ");
//			mpz_print_str_device(&_x1[0]);
//			printf(" ");
//			mpz_print_str_device(&_x2[0]);
//			printf("\n");
//		}
//		if(j == 1){
//			printf("mes2: ");
//			mpz_print_str_device(&_x1[0]);
//			printf(" ");
//			mpz_print_str_device(&_x2[0]);
//			printf("\n");
//		}

		for(int i = eLength - 2; i >= 0; i--){

			if(eBits[i] == 0){
				//x2 = _x1 * _x2
				mpz_mult(&tmp2[0], &_x1[0], &_x2[0]);
				//_x2 = REDC(rmod,n,n_,_x2,l)
				mpz_set( &_x2[0], REDC(rl, n, n_, &tmp2[0], &tmp[0], &t[0]) );
				//_x1 = _x1 * _x1
				mpz_set( &tmp[0], &_x1[0]);
				mpz_mult(&tmp2[0], &_x1[0], &tmp[0]);
				//_x1 = REDC(rmod,n,n_,_x1,l)
				mpz_set( &_x1[0], REDC(rl, n, n_, &tmp2[0], &tmp[0], &t[0]) );
			} else {
				//_x1 = _x1 * _x2
				mpz_mult(&tmp2[0], &_x1[0], &_x2[0]);
				//_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
				mpz_set( &_x1[0], REDC(rl, n, n_, &tmp2[0], &tmp[0], &t[0]) );
				//_x2 = _x2 * _x2
				mpz_set( &tmp[0], &_x2[0]);
				mpz_mult(&tmp2[0], &_x2[0], &tmp[0]);
				//_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
				mpz_set( &_x2[0], REDC(rl, n, n_, &tmp2[0], &tmp[0], &t[0]) );
			}

//			if(j == 0){
//				printf("mes1: ");
//				mpz_print_str_device(&_x1[0]);
//				printf(" ");
//				mpz_print_str_device(&_x2[0]);
//				printf("\n");
//			}
//			if(j == 1){
//				printf("mes2: ");
//				mpz_print_str_device(&_x1[0]);
//				printf(" ");
//				mpz_print_str_device(&_x2[0]);
//				printf("\n");
//			}
		}

		//_x1 = REDC(rmod,n,n_,_x1,l)
		mpz_set( &_x1[0], REDC(rl, n, n_, &_x1[0], &tmp[0], &t[0]) );

//		if(j == 0){
//			printf("mes1: ");
//			mpz_print_str_device(&_x1[0]);
//			printf("\n");
//		}
//		if(j == 1){
//			printf("mes2: ");
//			mpz_print_str_device(&_x1[0]);
//			printf("\n");
//		}
	}
}

void init(mpz_t* _x1, mpz_t* _x2, mpz_t* tmp, mpz_t* tmp2, mpz_t* t){

	mpz_init(&tmp[0]);////initial value not used
	mpz_init(&tmp2[0]);////initial value not used
	mpz_init(&_x1[0]);////initial value required
	mpz_init(&_x2[0]);////initial value required
	mpz_init(&t[0]);////initial value not used
}











