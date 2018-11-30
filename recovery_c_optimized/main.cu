#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.cu"
#include <time.h>
#include "cuda_mpz.h"
#include <gmp.h>

int CheckREDC(cuda_mpz_t* N, cuda_mpz_t* N_, cuda_mpz_t* T, cuda_mpz_t* tmp, cuda_mpz_t* t){

	//m = ((T & R) * N_) & R
	cuda_mpz_bitwise_truncate(t, T);
	cuda_mpz_mult(tmp, N_, t);
	cuda_mpz_bitwise_truncate_eq(tmp);

	//t = (T + m*N) >> L
	cuda_mpz_mult(t, tmp , N);
	cuda_mpz_add(tmp, T, t);
	cuda_mpz_bitwise_rshift(t, tmp);

	if (cuda_mpz_gte(t , N)){
		return 1;
    }
	else{
	    return 0;
	}
}

int CheckDivExp(cuda_mpz_t * mes1, cuda_mpz_t * mes2, int* eBits, int eLength, cuda_mpz_t* _x1_1, cuda_mpz_t* _x1_2, cuda_mpz_t* _x2_1, cuda_mpz_t* _x2_2,
		cuda_mpz_t* _x1_1_temp, cuda_mpz_t* _x1_2_temp, cuda_mpz_t* _x2_1_temp, cuda_mpz_t* _x2_2_temp,
		cuda_mpz_t* tmp_1, cuda_mpz_t* tmp_2, cuda_mpz_t* tmp2_1, cuda_mpz_t* tmp2_2, cuda_mpz_t* r2, cuda_mpz_t* n, cuda_mpz_t* n_,  cuda_mpz_t* t_1, cuda_mpz_t* t_2){

	//mes1 * r2
	cuda_mpz_mult(tmp2_1, mes1, r2);
	//mes2 * r2
	cuda_mpz_mult(tmp2_2, mes2, r2);
	//s1_1 = CheckREDC(rmod, n, n_, mes1 * r2, l)
	//s1_2 = CheckREDC(rmod, n, n_, mes2 * r2, l)
	int s1_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
	int s1_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

	//_x1_1 = REDC(rmod, n, n_, mes1 * r2, l)
	//_x1_2 = REDC(rmod, n, n_, mes2 * r2, l)
	cuda_mpz_set( _x1_1, REDC( n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x1_2, REDC( n, n_, tmp2_2, tmp_2, t_2) );
	//_x2_1 = _x1_1 * _x1_1
	//_x2_2 = _x1_2 * _x1_2
	cuda_mpz_mult(tmp2_1, _x1_1, t_1);
	cuda_mpz_mult(tmp2_2, _x1_2, t_2);

	//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int s2_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
	int s2_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);
	//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
	//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
	cuda_mpz_set( _x2_1, REDC( n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x2_2, REDC( n, n_, tmp2_2, tmp_2, t_2) );

	//for i in e_b[1:]:
	for(int i = 1; i < eLength; i++){ //big endian

		if(eBits[i] == 0){
			//_x2_1 = _x1_1 * _x2_1
			//_x2_2 = _x1_2 * _x2_2
			cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
			cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);

			//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			s2_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
			s2_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC( n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC( n, n_, tmp2_2, tmp_2, t_2) );

			//_x1_1 = _x1_1 * _x1_1
			//_x1_2 = _x1_2 * _x1_2
			cuda_mpz_set( tmp_1, _x1_1);
			cuda_mpz_mult(tmp2_1, _x1_1, tmp_1);
			cuda_mpz_set( tmp_2, _x1_2);
			cuda_mpz_mult(tmp2_2, _x1_2, tmp_2);

			//s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			//s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			s1_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
			s1_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

			//_x1_1 = REDC(rmod, n, n_, _x1_1, l)
			//_x1_2 = REDC(rmod, n, n_, _x1_2, l)
			cuda_mpz_set( _x1_1, REDC( n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x1_2, REDC( n, n_, tmp2_2, tmp_2, t_2) );
		} else{
			//_x1_1 = _x1_1 * _x2_1
			//_x1_2 = _x1_2 * _x2_2
			cuda_mpz_mult(tmp2_1, _x1_1, _x2_1);
			cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);

			//s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			//s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			s1_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
			s1_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

			//_x1_1 = REDC(rmod, n, n_, _x1_1, l)
			//_x1_2 = REDC(rmod, n, n_, _x1_2, l)
			cuda_mpz_set( _x1_1, REDC( n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x1_2, REDC( n, n_, tmp2_2, tmp_2, t_2) );

			//_x2_1 = _x2_1 * _x2_1
			//_x2_2 = _x2_2 * _x2_2
			cuda_mpz_set( tmp_1, _x2_1);
			cuda_mpz_mult(tmp2_1, _x2_1, tmp_1);
			cuda_mpz_set( tmp_2, _x2_2);
			cuda_mpz_mult(tmp2_2, _x2_2, tmp_2);

			//s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			//s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			s2_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
			s2_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC( n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC( n, n_, tmp2_2, tmp_2, t_2) );
		}
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
	int d0_s2_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
	//_x1_1 = _x1_1 * _x1_1
	cuda_mpz_set( tmp_1, _x1_1);
	cuda_mpz_mult(tmp2_1, _x1_1, tmp_1);
	//d0_s1_1 = CheckREDC(rmod, n, n_, _x1_1 ,l)
	int d0_s1_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);

	//_x2_2 = _x1_2 * _x2_2
	cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);
	//d0_s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int d0_s2_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);
	//_x1_2 = _x1_2 * _x1_2
	cuda_mpz_set( tmp_2, _x1_2);
	cuda_mpz_mult(tmp2_2, _x1_2, tmp_2);
	//d0_s1_2 = CheckREDC(rmod, n, n_, _x1_2 ,l)
	int d0_s1_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

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
	int d1_s1_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);
	//_x2_1 = _x2_1 * _x2_1
	cuda_mpz_set( tmp_1, _x2_1);
	cuda_mpz_mult(tmp2_1, _x2_1, tmp_1);
	//d1_s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	int d1_s2_1 = CheckREDC( n, n_, tmp2_1, tmp_1, t_1);

	//_x1_2 = _x1_2 * _x2_2
	cuda_mpz_mult(tmp2_2, _x1_2, _x2_2);
	//d1_s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
	int d1_s1_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);
	//_x2_2 = _x2_2 * _x2_2
	cuda_mpz_set( tmp_2, _x2_2);
	cuda_mpz_mult(tmp2_2, _x2_2, tmp_2);
	//d1_s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	int d1_s2_2 = CheckREDC( n, n_, tmp2_2, tmp_2, t_2);

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

	///////input control
	if (argc < 2){
		printf("sample size required.\n");
		exit(EXIT_FAILURE);
	}

	int peak_clk = 1;//kHz
	int dev_id = 0;
	cudaDeviceGetAttribute(&peak_clk, cudaDevAttrClockRate, dev_id);
	float clock_rate = (float) peak_clk;
	printf("clock_rate_out_kernel:%f\n", clock_rate);

	long x = strtol(argv[1], NULL, 10);
	long long unsigned pairs = x;
	unsigned thread_num = 2;
	long long unsigned data_num = pairs * thread_num;

	///////host memory
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc( 2 * sizeof(long long int));

	cuda_mpz_t h_n;
	cuda_mpz_t h_n_;
	cuda_mpz_t h_r2;

	///////get n
	char n_input[] = "00000038f6e8cfba55dd0e47";
	cuda_mpz_set_str_host(&h_n, n_input);

	///////get n_
	char n__input[] = "0000002e8457440e0d93c489";
	cuda_mpz_set_str_host(&h_n_, n__input);

	///////get r2
	char r2_input[] = "0000003709d17d8f8686609f";
	cuda_mpz_set_str_host(&h_r2, r2_input);

	///////get d
	char d_input[] = "1011011001001001010011110110010101010111001010110101111000111100001";
	int d_bitsLength = (int)strlen(d_input);
	int* dBits = (int *) malloc(sizeof(int) * d_bitsLength);
	int* dBits_d;
	cudaMalloc((void **) &dBits_d, sizeof(int) * d_bitsLength);

	int d_iterator = 0;
	while ( d_iterator < d_bitsLength){
		if( d_input[d_iterator] == '1'){//big endian
			dBits[d_iterator] = 1;
		}
		else{
			dBits[d_iterator] = 0;
		}
		d_iterator++;
	}
	cudaMemcpy(dBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);

	long long int *clockTable_d;
	cudaMalloc((void **) &clockTable_d, 2 * sizeof(long long int));

	///////get Messages
	long long unsigned mesSize = sizeof(cuda_mpz_t) * data_num;
	cuda_mpz_t *myMes1_h;
	myMes1_h = (cuda_mpz_t*) malloc (mesSize * 2); //CPU, bit1_div and bit0_div lists

	cuda_mpz_t *myMes1_d;
	cudaMalloc((cuda_mpz_t **) &myMes1_d, mesSize * 2); //GPU

	///////gen_pairs variables
	int	bit1_div_num = 0;
	int	bit0_div_num = 0;

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

	int known_bits[2048];
	known_bits[0] = 1;//first bit is always 1
	int known_bits_length = 1;
	int div_con = 0;
	int wrong_key = 0;

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

	printf("current bits: ");
	for(int i = 0; i < known_bits_length; i++){
		printf("%d", known_bits[i]);
	}
	printf("\n");

	while(known_bits_length < d_bitsLength - 1){

		bit1_div_num = 0;
		bit0_div_num = 0;

		while(1){
			mpz_urandomm (rand_num, rand_state, mod);
			cuda_mpz_set_gmp(&r1, rand_num);
			mpz_urandomm (rand_num, rand_state, mod);
			cuda_mpz_set_gmp(&r2, rand_num);

			div_con = CheckDivExp(&r1, &r2, known_bits, known_bits_length, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
											&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
											&tmp_1, &tmp_2, &tmp2_1, &tmp2_2,  &h_r2, &h_n, &h_n_,  &t_1, &t_2);

			if (div_con == 1 && bit1_div_num < data_num){
				cuda_mpz_set( &myMes1_h[bit1_div_num], &r1);
				bit1_div_num++;
				cuda_mpz_set( &myMes1_h[bit1_div_num], &r2);
				bit1_div_num++;
			}
			if (div_con == 4 && bit0_div_num < data_num){
				cuda_mpz_set( &myMes1_h[bit0_div_num + data_num], &r1);
				bit0_div_num++;
				cuda_mpz_set( &myMes1_h[bit0_div_num + data_num], &r2);
				bit0_div_num++;
			}
			if (bit1_div_num == data_num && bit0_div_num == data_num){
				break;
			}
		}

		long long int sum1 = 0;
		long long int sum4 = 0;

		////////////////////////////////////////////////////////////////converge for bit 0, diverge for bit 1
		cudaMemcpy(myMes1_d, myMes1_h, mesSize * 4 , cudaMemcpyHostToDevice);

		struct timespec ts1;/////////////////////////////////time
		clock_gettime(CLOCK_REALTIME, &ts1);/////////////////////////////////time

		printf("debug0\n");

		MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d);/////////////////////////////////////////kernel
		cudaDeviceSynchronize();

		printf("debug1\n");

		struct timespec ts2;/////////////////////////////////time
		clock_gettime(CLOCK_REALTIME, &ts2);/////////////////////////////////time
		long long unsigned time_interval = time_diff(ts1, ts2);/////////////////////////////////time
		printf("%lluns %fms %fs\n", time_interval,  ((double) time_interval) / 1000000,  ((double) time_interval) / 1000000000);/////////////////////////////////time

		cudaMemcpy(clockTable_h, clockTable_d, 2 * sizeof(long long int), cudaMemcpyDeviceToHost);

		sum1 = clockTable_h[0];
		sum1 = sum1 / pairs;
		sum4 = clockTable_h[1] - clockTable_h[0];
		sum4 = sum4 / pairs;

		long long int diff3 = sum1 - sum4;

		printf("%lld %lld %lld\n", sum1, sum4, diff3);
		printf ("bit1_div: %fms %lldcycles\n", sum1 / (float)clock_rate, sum1);
		printf ("bit0_div: %fms %lldcycles\n", sum4 / (float)clock_rate, sum4);

		if(diff3 > 1000){//bit is 1
			known_bits[known_bits_length] = 1;
			printf("bit is 1.\n");
		}else if(diff3 < -1000){//bit is 0
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

		if(known_bits[known_bits_length - 1] != dBits[known_bits_length - 1]){
			wrong_key = 1;
			printf("wrong key!");
			break;
		}
	}

	if(wrong_key == 0){
		known_bits[known_bits_length] = 1;//last bit is always 1
		printf("bit is 1.\n");

		known_bits_length++;

		printf("current bits: ");
		for(int i = 0; i < known_bits_length; i++){
			printf("%d", known_bits[i]);
		}
		printf("\n");
	}

	///////gmp clear
	gmp_randclear(rand_state);
	mpz_clear(rand_num);
	mpz_clear(mod);

	////////free device
	cudaFree(clockTable_d);
	cudaFree(dBits_d);
	cudaFree(myMes1_d);

	////////free host
	free(clockTable_h);
	free(myMes1_h);
	free(dBits);

    return 0;
}

