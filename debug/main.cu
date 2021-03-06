#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.cu"
#include <time.h>
#include "cuda_mpz.h"
#include <gmp.h>

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

int CheckDivExp(cuda_mpz_t * mes1, cuda_mpz_t * mes2, int* eBits, int eLength, cuda_mpz_t* _x1_1, cuda_mpz_t* _x1_2, cuda_mpz_t* _x2_1, cuda_mpz_t* _x2_2,
		cuda_mpz_t* _x1_1_temp, cuda_mpz_t* _x1_2_temp, cuda_mpz_t* _x2_1_temp, cuda_mpz_t* _x2_2_temp,
		cuda_mpz_t* tmp_1, cuda_mpz_t* tmp_2, cuda_mpz_t* tmp2_1, cuda_mpz_t* tmp2_2, int rl, cuda_mpz_t* r2, cuda_mpz_t* n, cuda_mpz_t* n_,  cuda_mpz_t* t_1, cuda_mpz_t* t_2){

	//mes1 * r2
	cuda_mpz_mult(tmp2_1, mes1, r2);
	//mes2 * r2
	cuda_mpz_mult(tmp2_2, mes2, r2);

	//s1_1 = CheckREDC(rmod, n, n_, mes1 * r2, l)
	//s1_2 = CheckREDC(rmod, n, n_, mes2 * r2, l)
	int s1_1 = CheckREDC(rl, n, n_, tmp2_1, tmp_1, t_1);
	int s1_2 = CheckREDC(rl, n, n_, tmp2_2, tmp_2, t_2);

	if (s1_1 != s1_2){ //previous bits are all convergent
		return 0;
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

	if (s2_1 != s2_2){ //previous bits are all convergent
		return 0;
	}

	//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
	//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
	cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
	cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

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

			if (s2_1 != s2_2){ //previous bits are all convergent
				return 0;
			}

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );

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

			if (s1_1 != s1_2){ //previous bits are all convergent
				return 0;
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

			if (s1_1 != s1_2){ //previous bits are all convergent
				return 0;
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

			if (s2_1 != s2_2){ //previous bits are all convergent
				return 0;
			}

			//_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			//_x2_2 = REDC(rmod, n, n_, _x2_2, l)
			cuda_mpz_set( _x2_1, REDC(rl, n, n_, tmp2_1, tmp_1, t_1) );
			cuda_mpz_set( _x2_2, REDC(rl, n, n_, tmp2_2, tmp_2, t_2) );
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

	if( d0_s1_1 != d0_s1_2 || d0_s2_1 != d0_s2_2 ){ //diverge for bit 0
		if ( d1_s1_1 != d1_s1_2 || d1_s2_1 != d1_s2_2 ){ //diverge for bit 0 and diverge for bit 1
			return 0
		} else{ //diverge for bit 0 and converge for bit 1
			return 4;
		}
	} else if ( d1_s1_1 != d1_s1_2 || d1_s2_1 != d1_s2_2 ){ //converge for bit 0, diverge for bit 1
		return 1
	}
	else{ //converge for bit 0 and converge for bit 1
		return 2
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

	struct timespec ts1;
	clock_gettime(CLOCK_REALTIME, &ts1);

	///////input control
	if (argc < 5){
		exit(EXIT_FAILURE);
	}

	long x = strtol(argv[5], NULL, 10);
	long long unsigned pairs = x;
	unsigned thread_num = 2;
	long long unsigned data_num = pairs * thread_num;

	///////host memory
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc( pairs * sizeof(long long int));

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
	cudaMemcpy(dBits_d, dBits, sizeof(int) * d_bitsLength, cudaMemcpyHostToDevice);

	///////device memory
	unsigned varSize = sizeof(cuda_mpz_t) * thread_num;

	long long int *clockTable_d;
	cuda_mpz_t *tmp;
	cuda_mpz_t *tmp2;
	cuda_mpz_t *d_t;
	cuda_mpz_t *_x1_cuda_mpz;
	cuda_mpz_t *_x2_cuda_mpz;
	cudaMalloc((void **) &clockTable_d, pairs * sizeof(long long int));
	cudaMalloc((void **) &tmp, varSize);
	cudaMalloc((void **) &tmp2, varSize);
	cudaMalloc((void **) &d_t, varSize);
	cudaMalloc((void **) &_x2_cuda_mpz, varSize);
	cudaMalloc((void **) &_x1_cuda_mpz, varSize);

	init<<<1, thread_num>>>(_x1_cuda_mpz, _x2_cuda_mpz, tmp, tmp2, d_t);
	cudaDeviceSynchronize();

	///////get Messages
	long long unsigned mesSize = sizeof(cuda_mpz_t) * data_num;
	cuda_mpz_t *myMes1_h;
	myMes1_h = (cuda_mpz_t*) malloc (mesSize);//CPU

	for(long long unsigned i = 0; i < data_num; i++){
		cuda_mpz_init(&myMes1_h[i]);
	}

	cuda_mpz_t *myMes1_d;
	cudaMalloc((cuda_mpz_t **) &myMes1_d, mesSize);//GPU


	///////gen_pairs variables
	int	bit1_div_num = pairs;
	int nondiv_num = pairs;
	int	bit0_div_num = pairs;

	cuda_mpz_t r1, r2;
	cuda_mpz_t _x1_1, _x1_2, _x2_1, _x2_2
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
	int known_bits_length = 1;
	int div_con = 0;

	///////gmp init
	mpz_t mod;
	mpz_t rand_num;
	mpz_init (mod);
	mpz_init (rand_num);

	mpz_set_str (mod, n_input, 16);

	///////RNG init
	gmp_randstate_t rand_state;
	gmp_randinit_default (rand_state);
	gmp_randseed_ui (rand_state, time(NULL));


	while(1){

		mpz_urandomm (rand_num, rand_state, mod);
		cuda_mpz_set_gmp(&r1, rand_num);
		mpz_urandomm (rand_num, rand_state, mod);
		cuda_mpz_set_gmp(&r2, rand_num);

		div_con = CheckDivExp(&r1, &r2, known_bits, known_bits_length, &_x1_1, &_x1_2, &_x2_1, &_x2_2,
				&_x1_1_temp, &_x1_2_temp, &_x2_1_temp, &_x2_2_temp,
				&tmp_1, &tmp_2, &tmp2_1, &tmp2_2, rl, &h_r2, &h_n, &h_n_,  &t_1, &t_2);

		if (div_con == 1 && bit1_div_num > 0){
			f1.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			bit1_div_num--;
		}
		if (div_con == 2 && nondiv_num > 0){
			f2.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			nondiv_num--;
		}
		if (div_con == 4 && bit0_div_num > 0){
			f3.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			bit0_div_num--;
		}
		if (bit1_div_num == 0 && nondiv_num == 0 && bit0_div_num == 0){
			break;
		}
	}

	///////gmp clear
	gmp_randclear (rand_state);
	mpz_clear (rand_num1);
	mpz_clear (mod);

	///////get Message pairs
	char* line = NULL;
	size_t len = 0;

	for(int i = 2; i < 5; i++){
		FILE* fp = fopen(argv[i], "r");//input from pair storage
		if (fp == NULL){
			exit(EXIT_FAILURE);
		}

		long long unsigned line_num = 0;
		while ((getline(&line, &len, fp)) != -1) {
			line[strcspn(line, "\n")] = 0;
			cuda_mpz_set_str_host(&myMes1_h[line_num], line);
			line_num++;
			if(line_num == data_num){
				break;
			}
		}
		fclose(fp);

		cudaMemcpy(myMes1_d, myMes1_h, mesSize, cudaMemcpyHostToDevice);

		MontSQMLadder<<<1, thread_num>>>(myMes1_d, pairs, _x1_cuda_mpz, _x2_cuda_mpz, tmp, tmp2, rl, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d, d_t);/////////////////////////////////////////kernel
		cudaDeviceSynchronize();

		cudaMemcpy(clockTable_h, clockTable_d, pairs * sizeof(long long int), cudaMemcpyDeviceToHost);

		long long unsigned sum = 0;

		for (long long unsigned q = 0; q < pairs; q++){
			sum += clockTable_h[q];
		}

		sum = sum / pairs;
		printf("%llu ", sum);
	}

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
	if (line){
		free(line);
	}
	free(clockTable_h);
	free(myMes1_h);
	free(eBits);
	free(dBits);

	/////////////////////////////////time
	struct timespec ts2;
	clock_gettime(CLOCK_REALTIME, &ts2);

	printf("%llu ", time_diff(ts1, ts2));

    return 0;
}

