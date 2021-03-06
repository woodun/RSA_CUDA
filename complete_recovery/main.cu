#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.cu"
#include <time.h>
#include "cuda_mpz.h"
#include <gmp.h>

//nvprof --print-gpu-trace --log-file prof1.txt ./main 1000 1 > nvprof.txt
int CheckREDC(int RL, cuda_mpz_t* N, cuda_mpz_t* N_, cuda_mpz_t* T, cuda_mpz_t* tmp, cuda_mpz_t* t){

	//m = ((T & R) * N_) & R
	cuda_mpz_bitwise_truncate(t, T, RL);
	cuda_mpz_mult(tmp, N_, t);
	cuda_mpz_bitwise_truncate_eq(tmp, RL);

	//t = (T + m*N) >> L
	cuda_mpz_mult(t, tmp , N);
	cuda_mpz_add(tmp, T, t);
	cuda_mpz_bitwise_rshift(t, tmp, RL);

	if (cuda_mpz_gte(t , N)){
		return 1;
    }
	else{
	    return 0;
	}
}

int Exp(cuda_mpz_t * mes1, cuda_mpz_t * mes2, int* eBits, int eLength, cuda_mpz_t* _x1_1, cuda_mpz_t* _x1_2, cuda_mpz_t* _x2_1, cuda_mpz_t* _x2_2,
		cuda_mpz_t* _x1_1_temp, cuda_mpz_t* _x1_2_temp, cuda_mpz_t* _x2_1_temp, cuda_mpz_t* _x2_2_temp,
		cuda_mpz_t* tmp_1, cuda_mpz_t* tmp_2, cuda_mpz_t* tmp2_1, cuda_mpz_t* tmp2_2, int rl, cuda_mpz_t* r2, cuda_mpz_t* n, cuda_mpz_t* n_,  cuda_mpz_t* t_1, cuda_mpz_t* t_2, long check_pre){

	int div_count = 0;

//	char test_str[1024];
//	printf("mes1: ");
//	printf("%s\n", cuda_mpz_get_str(mes1, test_str, 1024));
//	printf("mes2: ");
//	printf("%s\n", cuda_mpz_get_str(mes2, test_str, 1024));

	//mes1 * r2
	cuda_mpz_mult(tmp2_1, mes1, r2);
	//mes2 * r2
	cuda_mpz_mult(tmp2_2, mes2, r2);

	//s1_1 = CheckREDC(rmod, n, n_, mes1 * r2, l)
	//s1_2 = CheckREDC(rmod, n, n_, mes2 * r2, l)
	int s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	int s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	if (s1_1 != s1_2){
		div_count++;
	}

	//_x1_1 = REDC(rmod, n, n_, mes1 * r2, l)
	//_x1_2 = REDC(rmod, n, n_, mes2 * r2, l)
	cuda_mpz_set( _x1_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x1_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

	//_x2_1 = _x1_1 * _x1_1
	//_x2_2 = _x1_2 * _x1_2
	cuda_mpz_mult(tmp2_1, _x1_1, t_1);
	cuda_mpz_mult(tmp2_2, _x1_2, t_2);

	//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	int s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	if (s2_1 != s2_2){
		div_count++;
	}

	//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
	//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
	cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

//	printf("mes1: ");
//	printf("%s", cuda_mpz_get_str(_x1_1, test_str, 1024));
//	printf(" %s\n", cuda_mpz_get_str(_x2_1, test_str, 1024));
//	printf("mes2: ");
//	printf("%s", cuda_mpz_get_str(_x1_2, test_str, 1024));
//	printf(" %s\n", cuda_mpz_get_str(_x2_2, test_str, 1024));
//	printf("%d %d %d %d\n", s1_1, s1_2, s2_1, s2_2);

	//for i in e_b[1:]:
	for(int i = 1; i < eLength; i++){ //big endian

//		printf("ebit: %d\n", eBits[i]);
		if(eBits[i] == 0){
			//_x2_1 = _x1_1 * _x2_1
			//_x2_2 = _x1_2 * _x2_2
			cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
			cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);

			//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s2_1 != s2_2){
				div_count++;
			}

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

			//_x1_1 = _x1_1 * _x1_1
			//_x1_2 = _x1_2 * _x1_2
			cuda_mpz_set( tmp_1, _x1_1);
			cuda_mpz_mult(tmp2_1, _x1_1, tmp_1);
			cuda_mpz_set( tmp_2, _x1_2);
			cuda_mpz_mult(tmp2_2, _x1_2, tmp_2);

			//s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			//s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s1_1 != s1_2){
				div_count++;
			}

			//_x1_1 = REDC(rmod, n, n_, _x1_1, l)
			//_x1_2 = REDC(rmod, n, n_, _x1_2, l)
			cuda_mpz_set( _x1_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x1_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );
		} else{
			//_x1_1 = _x1_1 * _x2_1
			//_x1_2 = _x1_2 * _x2_2
			cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
			cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);

			//s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			//s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s1_1 != s1_2){
				div_count++;
			}

			//_x1_1 = REDC(rmod, n, n_, _x1_1, l)
			//_x1_2 = REDC(rmod, n, n_, _x1_2, l)
			cuda_mpz_set( _x1_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x1_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

			//_x2_1 = _x2_1 * _x2_1
			//_x2_2 = _x2_2 * _x2_2
			cuda_mpz_set( tmp_1, _x2_1);
			cuda_mpz_mult(tmp2_1, _x2_1, tmp_1);
			cuda_mpz_set( tmp_2, _x2_2);
			cuda_mpz_mult(tmp2_2, _x2_2, tmp_2);


			//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s2_1 != s2_2){
				div_count++;
			}

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );
		}

//		printf("mes1: ");
//		printf("%s", cuda_mpz_get_str(_x1_1, test_str, 1024));
//		printf(" %s\n", cuda_mpz_get_str(_x2_1, test_str, 1024));
//		printf("mes2: ");
//		printf("%s", cuda_mpz_get_str(_x1_2, test_str, 1024));
//		printf(" %s\n", cuda_mpz_get_str(_x2_2, test_str, 1024));
//
//		printf("%d %d %d %d\n", s1_1, s1_2, s2_1, s2_2);
	}

    //s1 = CheckREDC(rmod,n,n_,_x1_1,l)
    //s2 = CheckREDC(rmod,n,n_,_x1_2,l)
	int s1 = CheckREDC(rl, n, n_, _x1_1, tmp_1, t_1);
	int s2 = CheckREDC(rl, n, n_, _x1_2, tmp_2, t_2);

    //_x1_1 = REDC(rmod,n,n_,_x1_1,l)
    //_x1_2 = REDC(rmod,n,n_,_x1_2,l)
	cuda_mpz_set( _x1_1, REDC(rl, n, n_, _x1_1, tmp_1, t_1) );
	cuda_mpz_set( _x1_2, REDC(rl, n, n_, _x1_2, tmp_2, t_2) );

//	printf("mes1: ");
//	printf("%s\n", cuda_mpz_get_str(_x1_1, test_str, 1024));
//	printf("mes2: ");
//	printf("%s\n", cuda_mpz_get_str(_x1_2, test_str, 1024));
//
//	printf("%d %d\n", s1, s2);

	if (s1 != s2){
		div_count++;
	}

//	printf("%d\n", div_count);

	return div_count;
}

int CheckDivExp(cuda_mpz_t * mes1, cuda_mpz_t * mes2, int* eBits, int eLength, cuda_mpz_t* _x1_1, cuda_mpz_t* _x1_2, cuda_mpz_t* _x2_1, cuda_mpz_t* _x2_2,
		cuda_mpz_t* _x1_1_temp, cuda_mpz_t* _x1_2_temp, cuda_mpz_t* _x2_1_temp, cuda_mpz_t* _x2_2_temp,
		cuda_mpz_t* tmp_1, cuda_mpz_t* tmp_2, cuda_mpz_t* tmp2_1, cuda_mpz_t* tmp2_2, int rl, cuda_mpz_t* r2, cuda_mpz_t* n, cuda_mpz_t* n_,  cuda_mpz_t* t_1, cuda_mpz_t* t_2, long check_pre){

	int div_count = 0;

//	char test_str[1024];
//	printf("mes1: ");
//	printf("%s\n", cuda_mpz_get_str(mes1, test_str, 1024));
//	printf("mes2: ");
//	printf("%s\n", cuda_mpz_get_str(mes2, test_str, 1024));

	//mes1 * r2
	cuda_mpz_mult(tmp2_1, mes1, r2);
	//mes2 * r2
	cuda_mpz_mult(tmp2_2, mes2, r2);

	//s1_1 = CheckREDC(rmod, n, n_, mes1 * r2, l)
	//s1_2 = CheckREDC(rmod, n, n_, mes2 * r2, l)
	int s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	int s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	if (s1_1 != s1_2){
//		printf("asd\n");
//		return 0;
		div_count++;
	}

	//_x1_1 = REDC(rmod, n, n_, mes1 * r2, l)
	//_x1_2 = REDC(rmod, n, n_, mes2 * r2, l)
	cuda_mpz_set( _x1_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x1_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

	//_x2_1 = _x1_1 * _x1_1
	//_x2_2 = _x1_2 * _x1_2
	cuda_mpz_mult(tmp2_1, _x1_1, t_1);
	cuda_mpz_mult(tmp2_2, _x1_2, t_2);

	//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	int s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	if (s2_1 != s2_2){
//		printf("sdf\n");
//		return 0;
		div_count++;
	}

	//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
	//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
	cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

//	printf("mes1: ");
//	printf("%s", cuda_mpz_get_str(_x1_1, test_str, 1024));
//	printf(" %s\n", cuda_mpz_get_str(_x2_1, test_str, 1024));
//	printf("mes2: ");
//	printf("%s", cuda_mpz_get_str(_x1_2, test_str, 1024));
//	printf(" %s\n", cuda_mpz_get_str(_x2_2, test_str, 1024));

	//for i in e_b[1:]:
	for(int i = 1; i < eLength; i++){ //big endian

		if(eBits[i] == 0){
			//_x2_1 = _x1_1 * _x2_1
			//_x2_2 = _x1_2 * _x2_2
			cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
			cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);

			//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s2_1 != s2_2){
//				return 0;
				div_count++;
			}

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

			//_x1_1 = _x1_1 * _x1_1
			//_x1_2 = _x1_2 * _x1_2
			cuda_mpz_set( tmp_1, _x1_1);
			cuda_mpz_mult(tmp2_1, _x1_1, tmp_1);
			cuda_mpz_set( tmp_2, _x1_2);
			cuda_mpz_mult(tmp2_2, _x1_2, tmp_2);

			//s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			//s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s1_1 != s1_2){
//				return 0;
				div_count++;
			}

			//_x1_1 = REDC(rmod, n, n_, _x1_1, l)
			//_x1_2 = REDC(rmod, n, n_, _x1_2, l)
			cuda_mpz_set( _x1_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x1_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );
		} else{
			//_x1_1 = _x1_1 * _x2_1
			//_x1_2 = _x1_2 * _x2_2
			cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
			cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);

			//s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			//s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s1_1 != s1_2){
//				return 0;
				div_count++;
			}

			//_x1_1 = REDC(rmod, n, n_, _x1_1, l)
			//_x1_2 = REDC(rmod, n, n_, _x1_2, l)
			cuda_mpz_set( _x1_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x1_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

			//_x2_1 = _x2_1 * _x2_1
			//_x2_2 = _x2_2 * _x2_2
			cuda_mpz_set( tmp_1, _x2_1);
			cuda_mpz_mult(tmp2_1, _x2_1, tmp_1);
			cuda_mpz_set( tmp_2, _x2_2);
			cuda_mpz_mult(tmp2_2, _x2_2, tmp_2);


			//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
			s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

			if (s2_1 != s2_2){
//				return 0;
				div_count++;
			}

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );
		}

//		printf("mes1: ");
//		printf("%s", cuda_mpz_get_str(_x1_1, test_str, 1024));
//		printf(" %s\n", cuda_mpz_get_str(_x2_1, test_str, 1024));
//		printf("mes2: ");
//		printf("%s", cuda_mpz_get_str(_x1_2, test_str, 1024));
//		printf(" %s\n", cuda_mpz_get_str(_x2_2, test_str, 1024));
	}

	if(div_count != eLength && check_pre == 1){ //total divergence number
		return 0;
	}

	//_x1_1_temp = _x1_1
	cuda_mpz_set( _x1_1_temp, _x1_1);
	//_x2_1_temp = _x2_1
	cuda_mpz_set( _x2_1_temp, _x2_1);
	//_x1_2_temp = _x1_2
	cuda_mpz_set( _x1_2_temp, _x1_2);
	//_x2_2_temp = _x2_2
	cuda_mpz_set( _x2_2_temp, _x2_2);

	//simulate exp bit 0
	//_x2_1 = _x1_1 * _x2_1
	cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
	//d0_s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	int d0_s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	//_x1_1 = _x1_1 * _x1_1
	cuda_mpz_set( tmp_1, _x1_1);
	cuda_mpz_mult(tmp2_1, _x1_1, tmp_1);
	//d0_s1_1 = CheckREDC(rmod, n, n_, _x1_1 ,l)
	int d0_s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);

	//_x2_2 = _x1_2 * _x2_2
	cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);
	//d0_s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int d0_s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);
	//_x1_2 = _x1_2 * _x1_2
	cuda_mpz_set( tmp_2, _x1_2);
	cuda_mpz_mult(tmp2_2, _x1_2, tmp_2);
	//d0_s1_2 = CheckREDC(rmod, n, n_, _x1_2 ,l)
	int d0_s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	//simulate exp bit 1
	//_x1_1 = _x1_1_temp
	cuda_mpz_set( _x1_1, _x1_1_temp);
	//_x2_1 = _x2_1_temp
	cuda_mpz_set( _x2_1, _x2_1_temp);
	//_x1_2 = _x1_2_temp
	cuda_mpz_set( _x1_2, _x1_2_temp);
	//_x2_2 = _x2_2_temp
	cuda_mpz_set( _x2_2, _x2_2_temp);

	//_x1_1 = _x1_1 * _x2_1
	cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
	//d1_s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
	int d1_s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	//_x2_1 = _x2_1 * _x2_1
	cuda_mpz_set( tmp_1, _x2_1);
	cuda_mpz_mult(tmp2_1, _x2_1, tmp_1);
	//d1_s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	int d1_s2_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);

	//_x1_2 = _x1_2 * _x2_2
	cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);
	//d1_s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
	int d1_s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);
	//_x2_2 = _x2_2 * _x2_2
	cuda_mpz_set( tmp_2, _x2_2);
	cuda_mpz_mult(tmp2_2, _x2_2, tmp_2);
	//d1_s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int d1_s2_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	if ( (d0_s1_1 != d0_s1_2 && d0_s2_1 == d0_s2_2) || (d0_s1_1 == d0_s1_2 && d0_s2_1 != d0_s2_2) ){ //diverge for bit 0 (1 0) or (0 1)
		if ( (d1_s1_1 != d1_s1_2 && d1_s2_1 == d1_s2_2) or (d1_s1_1 == d1_s1_2 && d1_s2_1 != d1_s2_2) ){ //diverge for bit 0, diverge for bit 1 (1 0) or (0 1)
			//printf ("debug3\n");
			return 3;
		} else if ( d1_s1_1 == d1_s1_2 && d1_s2_1 == d1_s2_2 ) { //diverge for bit 0, converge for bit 1 (0 0)
			//printf ("debug4\n");
			return 4;
		} else {
			return 0;
		}
	} else if (d0_s1_1 == d0_s1_2 && d0_s2_1 == d0_s2_2 ){ //converge for bit 0 (0 0)
		if ( (d1_s1_1 != d1_s1_2 && d1_s2_1 == d1_s2_2) || (d1_s1_1 == d1_s1_2 && d1_s2_1 != d1_s2_2) ){ //converge for bit 0, diverge for bit 1 (1 0) or (0 1)
			//printf ("debug1\n");
			return 1;
		} else if ( d1_s1_1 == d1_s1_2 && d1_s2_1 == d1_s2_2 ){ //converge for bit 0, converge for bit 1 (0 0)
			//printf ("debug2\n");
			return 2;
		} else {
			return 0;
		}
	} else {
		return 0;
	}
//	if( d0_s1_1 != d0_s1_2 || d0_s2_1 != d0_s2_2 ){ //diverge for bit 0
//		if ( d1_s1_1 != d1_s1_2 || d1_s2_1 != d1_s2_2 ){ //diverge for bit 0 and diverge for bit 1
////			printf("debug0\n");
//			return 0;
//		} else{ //diverge for bit 0 and converge for bit 1
////			printf("debug3\n");
//			return 3;
//		}
//	} else if ( d1_s1_1 != d1_s1_2 || d1_s2_1 != d1_s2_2 ){ //converge for bit 0, diverge for bit 1
////		printf("debug1\n");
//		return 1;
//	} else{ //converge for bit 0 and converge for bit 1
////		printf("debug2\n");
//		return 2;
//	}
}

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
	if (argc < 3){
		exit(EXIT_FAILURE);
	}

	long x = strtol(argv[1], NULL, 10);
	long long unsigned pairs = x;
	unsigned thread_num = 2;
	long long unsigned data_num = pairs * thread_num;

	///////host memory
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc( 4 * pairs * sizeof(long long int));

	cuda_mpz_t h_n;
	cuda_mpz_t h_n_;
	cuda_mpz_t h_r2;
	int rl = 70;

	cuda_mpz_init(&h_n);
	cuda_mpz_init(&h_n_);
	cuda_mpz_init(&h_r2);

	///////get n
	char n_input[] = "00000038f6e8cfba55dd0e47";
	cuda_mpz_set_str_host(&h_n, n_input);

	///////get n_
	char n__input[] = "0000002e8457440e0d93c489";
	cuda_mpz_set_str_host(&h_n_, n__input);

	///////get r2
	char r2_input[] = "0000003709d17d8f8686609f";
	cuda_mpz_set_str_host(&h_r2, r2_input);

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
	char d_input[] = "1011011001001001010011110110010101010111001010110101111000111100001"; //big endian
	//char d_input[] = "1000100010110110111110111000110000000001011000001000011010101101000101"; //big endian

	int d_bitsLength = (int)strlen(d_input);

	int* dBits = (int *) malloc(sizeof(int) * d_bitsLength);

	int* dBits_d;
	cudaMalloc((void **) &dBits_d, sizeof(int) * d_bitsLength);

	int d_iterator = d_bitsLength - 1;
	while ( d_iterator > 0){
        if( d_input[d_bitsLength - 1 - d_iterator] == '1'){//little endian
            dBits[d_iterator] = 1;
        }
        else{
            dBits[d_iterator] = 0;
        }
        d_iterator--;
	}
	dBits[d_iterator] = 1;

//	int d_iterator = 0;
//	while ( d_iterator < d_bitsLength){
//        if( d_input[d_iterator] == '1'){//big endian
//            dBits[d_iterator] = 1;
//        }
//        else{
//            dBits[d_iterator] = 0;
//        }
//        d_iterator++;
//	}

	cudaMemcpy(dBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);

	///////device memory
	unsigned varSize = sizeof(cuda_mpz_t) * thread_num;

	long long int *clockTable_d;
	cuda_mpz_t *tmp;
	cuda_mpz_t *tmp2;
	cuda_mpz_t *d_t;
	cuda_mpz_t *_x1_cuda_mpz;
	cuda_mpz_t *_x2_cuda_mpz;
	cudaMalloc((void **) &clockTable_d, 4 * pairs * sizeof(long long int));
	cudaMalloc((void **) &tmp, varSize);
	cudaMalloc((void **) &tmp2, varSize);
	cudaMalloc((void **) &d_t, varSize);
	cudaMalloc((void **) &_x2_cuda_mpz, varSize);
	cudaMalloc((void **) &_x1_cuda_mpz, varSize);

	////////////////////////////////////////////////////////////////initialize
	init<<<1, thread_num>>>(_x1_cuda_mpz, _x2_cuda_mpz, tmp, tmp2, d_t);
	cudaDeviceSynchronize();

	///////get Messages
	long long unsigned mesSize = sizeof(cuda_mpz_t) * data_num;
	cuda_mpz_t *myMes1_h;
	myMes1_h = (cuda_mpz_t*) malloc (mesSize * 4); //CPU list converge for bit 0, diverge for bit 1
	//cuda_mpz_t *myMes2_h;
	//myMes2_h = (cuda_mpz_t*) malloc (mesSize); //CPU list converge for bit 0 and converge for bit 1
	//cuda_mpz_t *myMes3_h;
	//myMes3_h = (cuda_mpz_t*) malloc (mesSize); //CPU list diverge for bit 0 and converge for bit 1

	for(long long unsigned i = 0; i < data_num * 4; i++){
		cuda_mpz_init(&myMes1_h[i]);
		//cuda_mpz_init(&myMes2_h[i]);
		//cuda_mpz_init(&myMes3_h[i]);
	}

	cuda_mpz_t *myMes1_d;
	cudaMalloc((cuda_mpz_t **) &myMes1_d, mesSize * 4); //GPU

	///////gen_pairs variables
	int	bit1_div_num = 0;
	int nondiv_num = 0;
	int bothdiv_num = 0;
	int	bit0_div_num = 0;

//	int bit1_div_sum = 0; //0 1
//    int nondiv_sum = 0; //0 0
//    int bothdiv_sum = 0; //1 1
//    int bit0_div_sum = 0; //1 0

	cuda_mpz_t r1, r2;
	cuda_mpz_t _x1_1, _x1_2, _x2_1, _x2_2;
	cuda_mpz_t _x1_1_temp, _x1_2_temp, _x2_1_temp, _x2_2_temp;
	cuda_mpz_t tmp_1, tmp_2, tmp2_1, tmp2_2, t_1, t_2;

	cuda_mpz_init(&r1);
	cuda_mpz_init(&r2);
	cuda_mpz_init(&_x1_1);
	cuda_mpz_init(&_x1_2);
	cuda_mpz_init(&_x2_1);
	cuda_mpz_init(&_x2_2);
	cuda_mpz_init(&_x1_1_temp);
	cuda_mpz_init(&_x1_2_temp);
	cuda_mpz_init(&_x2_1_temp);
	cuda_mpz_init(&_x2_2_temp);
	cuda_mpz_init(&tmp_1);
	cuda_mpz_init(&tmp_2);
	cuda_mpz_init(&tmp2_1);
	cuda_mpz_init(&tmp2_2);
	cuda_mpz_init(&t_1);
	cuda_mpz_init(&t_2);

	long check_pre = strtol(argv[2], NULL, 10);
	int known_bits[2048];
	known_bits[0] = 1;//first bit is always 1
	//known_bits[1] = 0;
	//known_bits[2] = 1;
	int known_bits_length = 1;
	int total_bits_length = 67;
	int div_con = 0;

	///////gmp init
	mpz_t mod;
	mpz_t rand_num;
	mpz_init (mod);
	mpz_init (rand_num);

	mpz_set_str (mod, n_input, 16);

	///////RNG init
	gmp_randstate_t rand_state;
	//gmp_randinit_default (rand_state);
	gmp_randinit_mt(rand_state);
	gmp_randseed_ui (rand_state, time(NULL));
	//gmp_randseed_ui (rand_state, 0);

	//printf("debug1\n");

	printf("current bits: ");
	for(int i = 0; i < known_bits_length; i++){
		printf("%d", known_bits[i]);
	}
	printf("\n");

	while(known_bits_length < total_bits_length - 1){

		bit1_div_num = 0;
		nondiv_num = 0;
		bothdiv_num = 0;
		bit0_div_num = 0;

		while(1){

			mpz_urandomm (rand_num, rand_state, mod);
			cuda_mpz_set_gmp(&r1, rand_num);
			mpz_urandomm (rand_num, rand_state, mod);
			cuda_mpz_set_gmp(&r2, rand_num);

	//		char r1_str[] = "00000000167a024204f7c1bd89";
	//		cuda_mpz_set_str_host(&r1, r1_str);
	//		char r2_str[] = "0000000088133287637ebdcdb";
	//		cuda_mpz_set_str_host(&r2, r2_str);
	//
	//		char test_str[1024];
	//		printf("%s\n", cuda_mpz_get_str(&r1, test_str, 1024));
	//		printf("%s\n", cuda_mpz_get_str(&r2, test_str, 1024));

			div_con = CheckDivExp(&r1, &r2, known_bits, known_bits_length, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
											&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
											&tmp_1, &tmp_2, &tmp2_1, &tmp2_2, rl, &h_r2, &h_n, &h_n_,  &t_1, &t_2, check_pre);

			if (div_con == 1 && bit1_div_num < data_num){
				cuda_mpz_set( &myMes1_h[bit1_div_num], &r1);
				bit1_div_num++;
				cuda_mpz_set( &myMes1_h[bit1_div_num], &r2);
				bit1_div_num++;

	//			bit1_div_sum+=Exp(&r1, &r2, dBits, d_bitsLength, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
	//										&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
	//										&tmp_1, &tmp_2, &tmp2_1, &tmp2_2, rl, &h_r2, &h_n, &h_n_,  &t_1, &t_2, check_pre);
			}
			if (div_con == 2 && nondiv_num < data_num){
				cuda_mpz_set( &myMes1_h[nondiv_num + data_num], &r1);
				nondiv_num++;
				cuda_mpz_set( &myMes1_h[nondiv_num + data_num], &r2);
				nondiv_num++;

	//			nondiv_sum+=Exp(&r1, &r2, dBits, d_bitsLength, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
	//										&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
	//										&tmp_1, &tmp_2, &tmp2_1, &tmp2_2, rl, &h_r2, &h_n, &h_n_,  &t_1, &t_2, check_pre);

	//			break;
			}
			if (div_con == 3 && bothdiv_num < data_num){
				cuda_mpz_set( &myMes1_h[bothdiv_num + data_num * 2], &r1);
				bothdiv_num++;
				cuda_mpz_set( &myMes1_h[bothdiv_num + data_num * 2], &r2);
				bothdiv_num++;

	//			bothdiv_sum+=Exp(&r1, &r2, dBits, d_bitsLength, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
	//										&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
	//										&tmp_1, &tmp_2, &tmp2_1, &tmp2_2, rl, &h_r2, &h_n, &h_n_,  &t_1, &t_2, check_pre);
			}
			if (div_con == 4 && bit0_div_num < data_num){
				cuda_mpz_set( &myMes1_h[bit0_div_num + data_num * 3], &r1);
				bit0_div_num++;
				cuda_mpz_set( &myMes1_h[bit0_div_num + data_num * 3], &r2);
				bit0_div_num++;

	//			bit0_div_sum+=Exp(&r1, &r2, dBits, d_bitsLength, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
	//										&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
	//										&tmp_1, &tmp_2, &tmp2_1, &tmp2_2, rl, &h_r2, &h_n, &h_n_,  &t_1, &t_2, check_pre);
			}
	//		if (div_con == 2 && nondiv_num > 0){
	//			nondiv_num--;
	//			cuda_mpz_set( &myMes2_h[nondiv_num], &r2);
	//			nondiv_num--;
	//			cuda_mpz_set( &myMes2_h[nondiv_num], &r1);
	//		}
	//		if (div_con == 3 && bit0_div_num > 0){
	//			bit0_div_num--;
	//			cuda_mpz_set( &myMes3_h[bit0_div_num], &r2);
	//			bit0_div_num--;
	//			cuda_mpz_set( &myMes3_h[bit0_div_num], &r1);
	//		}
			if (bit1_div_num == data_num && nondiv_num == data_num && bothdiv_num == data_num && bit0_div_num == data_num){

	//			double sum1 = (double) bit1_div_sum / pairs;
	//			double sum2 = (double) nondiv_sum / pairs;
	//			double sum3 = (double) bothdiv_sum / pairs;
	//			double sum4 = (double) bit0_div_sum / pairs;
	//
	//			printf("#########################################CPU output###########################################\n");
	//			printf("%f\n", sum1);
	//			printf("%f\n", sum2);
	//			printf("%f\n", sum3);
	//			printf("%f\n", sum4);
	//
	//			double diff1 = sum1 - sum2; // close to zero means 0, greater than zero means 1
	//			double diff2 = sum4 - sum2; // close to zero means 1, greater than zero means 0
	//			double mean1 = (diff1 + diff2) / 2;
	//			printf("%f %f %f\n",diff1,diff2,mean1);
	//
	//			double diff3 = sum1 - sum3; // close to zero means 1, smaller than zero means 0
	//			double diff4 = sum4 - sum3; // close to zero means 0, smaller than zero means 1
	//			double mean2 = (diff3 + diff4) / 2;
	//			printf("%f %f %f\n", diff3,diff4, mean2);
	//
	//			double val1 = mean1 + mean2;
	//			double sign1 = mean1 * mean2;
	//			printf("%f %f\n", val1, sign1);
	//
	//			double diff5 = sum1 - sum4; // greater means 1
	//			double diff6 = sum3 - sum2; // must be greater
	//			printf("%f %f\n",diff5,diff6);

				break;
			}
		}

	//	exit(0);

		long long int sum1 = 0;
		long long int sum2 = 0;
		long long int sum3 = 0;
		long long int sum4 = 0;

		//printf("debug2\n");

		////////////////////////////////////////////////////////////////converge for bit 0, diverge for bit 1
		cudaMemcpy(myMes1_d, myMes1_h, mesSize * 4 , cudaMemcpyHostToDevice);
		MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs * 4, _x1_cuda_mpz, _x2_cuda_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
		cudaDeviceSynchronize();
		cudaMemcpy(clockTable_h, clockTable_d, 4 * pairs * sizeof(long long int), cudaMemcpyDeviceToHost);

		//printf("debug3\n");

		for (long long unsigned q = 0; q < pairs * 1; q++){
			sum1 += clockTable_h[q];
		}
		sum1 = sum1 / pairs;
		for (long long unsigned q = pairs * 1; q < pairs * 2; q++){
			sum2 += clockTable_h[q];
		}
		sum2 = sum2 / pairs;
		for (long long unsigned q = pairs * 2; q < pairs * 3; q++){
			sum3 += clockTable_h[q];
		}
		sum3 = sum3 / pairs;
		for (long long unsigned q = pairs * 3; q < pairs * 4; q++){
			sum4 += clockTable_h[q];
		}
		sum4 = sum4 / pairs;

	//	////////////////////////////////////////////////////////////////converge for bit 0 and converge for bit 1
	//	cudaMemcpy(myMes1_d, myMes2_h, mesSize, cudaMemcpyHostToDevice);
	//	MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, _x1_cuda_mpz, _x2_cuda_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
	//	cudaDeviceSynchronize();
	//	cudaMemcpy(clockTable_h, clockTable_d, pairs * sizeof(long long int), cudaMemcpyDeviceToHost);
	//
	//	for (long long unsigned q = 0; q < pairs; q++){
	//		sum2 += clockTable_h[q];
	//	}
	//	sum2 = sum2 / pairs;
	//
	//	////////////////////////////////////////////////////////////////diverge for bit 0 and converge for bit 1
	//	cudaMemcpy(myMes1_d, myMes3_h, mesSize, cudaMemcpyHostToDevice);
	//	MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, _x1_cuda_mpz, _x2_cuda_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
	//	cudaDeviceSynchronize();
	//	cudaMemcpy(clockTable_h, clockTable_d, pairs * sizeof(long long int), cudaMemcpyDeviceToHost);
	//
	//	for (long long unsigned q = 0; q < pairs; q++){
	//		sum3 += clockTable_h[q];
	//	}
	//	sum3 = sum3 / pairs;

		long long int diff1 = abs(sum1 - sum2);
		long long int diff2 = abs(sum2 - sum4);
		long long int diff3 = sum1 - sum4;

		printf("%lld %lld %lld %lld %lld %lld %f %f\n", sum1, sum2, sum3, diff1, diff2, diff3, ((double) diff1) / diff2, ((double) diff2) / diff1);

	//	if(((double) diff1) / diff2 > 1.3){//bit is 1
	//		printf("bit is 1.\n");
	//	}else if(((double) diff2) / diff1 > 1.3){//bit is 0
	//		printf("bit is 0.\n");
	//	}else{//EOB
	//		printf("end of bits.\n");
	//	}

		if(diff3 > 5000){//bit is 1
			known_bits[known_bits_length] = 1;
			printf("bit is 1.\n");
		}else if(diff3 < -5000){//bit is 0

			known_bits[known_bits_length] = 0;
			printf("bit is 0.\n");
		}else{//EOB
			//printf("end of bits.\n");
			printf("bit not accepted.\n");
			continue;
		}

		known_bits_length++;

		printf("current bits: ");
		for(int i = 0; i < known_bits_length; i++){
			printf("%d", known_bits[i]);
		}
		printf("\n");
	}

	known_bits[known_bits_length] = 1;//last bit is always 1
	printf("bit is 1.\n");

	known_bits_length++;

	printf("current bits: ");
	for(int i = 0; i < known_bits_length; i++){
		printf("%d", known_bits[i]);
	}
	printf("\n");

	///////gmp clear
	gmp_randclear (rand_state);
	mpz_clear (rand_num);
	mpz_clear (mod);

	////////free device
	cudaFree(clockTable_d);
	cudaFree(eBits_d);
	cudaFree(dBits_d);
	cudaFree(myMes1_d);
	cudaFree(tmp);
	cudaFree(tmp2);
	cudaFree(d_t);
	cudaFree(_x1_cuda_mpz);
	cudaFree(_x2_cuda_mpz);

	////////free host
	free(clockTable_h);
	free(myMes1_h);
	//free(myMes2_h);
	//free(myMes3_h);
	free(eBits);
	free(dBits);

	/////////////////////////////////time
	struct timespec ts2;
	clock_gettime(CLOCK_REALTIME, &ts2);

	printf("%llu ", time_diff(ts1, ts2));

    return 0;
}

