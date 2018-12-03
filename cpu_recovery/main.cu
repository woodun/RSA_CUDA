#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cuda_mpz.h"
#include <gmp.h>

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

__device__ __host__ inline cuda_mpz_t* CUDA_REDC(cuda_mpz_t* N, cuda_mpz_t* N_, cuda_mpz_t* T, cuda_mpz_t* tmp, cuda_mpz_t* t){//cuda_mpz_t* RMOD, int L, cuda_mpz_t* N, cuda_mpz_t* N_ should not be changed.

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

void MontSQMLadder(cuda_mpz_t * mes1, long long unsigned pairs, cuda_mpz_t r2, cuda_mpz_t vn, cuda_mpz_t vn_, int* eBits, int eLength, long long int* clockTable) {

	cuda_mpz_t tmp; //capacity will become a problem for shared memory with large keys
	cuda_mpz_t tmp2;
	cuda_mpz_t t;
	cuda_mpz_t _x1;
	cuda_mpz_t _x2;

//	long long int &t1, &t2;
//	struct &timespec &ts1;/////////////////////////////////time
//	clock_gettime(CLOCK_REALTIME, &ts1);/////////////////////////////////time

	long long int sum1 = 0;

	///////////////////initialize
	cuda_mpz_init(&tmp);
	cuda_mpz_init(&tmp2);
	cuda_mpz_init(&t);
	cuda_mpz_init(&_x1);
	cuda_mpz_init(&_x2);

	cuda_mpz_t* n = &vn;
	cuda_mpz_t* n_ = &vn_;

	//to accelerate &the experiment, we put all messages in one kernel launch. In &the real case, each message causes one kernel launch.
	for(unsigned iter1 = 0; iter1 < 2 * pairs; ++iter1){

//		t1 = clock64();//beginning of necessary instructions within &the kernel

		cuda_mpz_set(&_x1, &mes1[iter1]);//next &_x1 access will cause L1 miss if &the L1 policy is write evict, same as using mutiple kernels.

		//&_x1 = CUDA_REDC(rmod,n,n_,mes*r2,l)
		cuda_mpz_mult( &tmp2, &_x1, &r2);
		cuda_mpz_set( &_x1, CUDA_REDC(n, n_, &tmp2, &tmp, &t) );
		//x2 = &_x1 * &_x1
		cuda_mpz_mult( &tmp2, &_x1, &t);
		//&_x2 = CUDA_REDC(rmod,n,n_,&_x2,l)
		cuda_mpz_set( &_x2, CUDA_REDC(n, n_, &tmp2, &tmp, &t) );

		for(int i = 1; i < eLength; ++i){
			if(eBits[i] == 0){
				//x2 = &_x1 * &_x2
				cuda_mpz_mult( &tmp2, &_x1,  &_x2);
				//&_x2 = CUDA_REDC(rmod,n,n_,&_x2,l)
				cuda_mpz_set(  &_x2, CUDA_REDC(n, n_, &tmp2, &tmp, &t) );
				//&_x1 = &_x1 * &_x1
				cuda_mpz_set( &tmp, &_x1);
				cuda_mpz_mult( &tmp2, &_x1, &tmp);
				//&_x1 = CUDA_REDC(rmod,n,n_,&_x1,l)
				cuda_mpz_set( &_x1, CUDA_REDC(n, n_, &tmp2, &tmp, &t) );
			} else {
				//&_x1 = &_x1 * &_x2
				cuda_mpz_mult( &tmp2, &_x1,  &_x2);
				//&_x1 = CUDA_REDC(rmod,n,n_,&_x1,l) #changes: more efficient
				cuda_mpz_set( &_x1, CUDA_REDC(n, n_, &tmp2, &tmp, &t) );
				//&_x2 = &_x2 * &_x2
				cuda_mpz_set( &tmp,  &_x2);
				cuda_mpz_mult( &tmp2,  &_x2, &tmp);
				//&_x2 = CUDA_REDC(rmod,n,n_,&_x2,l) #changes: more efficient
				cuda_mpz_set(  &_x2, CUDA_REDC(n, n_, &tmp2, &tmp, &t) );
			}
		}

		//&_x1 = CUDA_REDC(rmod,n,n_,&_x1,l)
		cuda_mpz_set( &_x1, CUDA_REDC(n, n_, &_x1, &tmp, &t) );

//		t2 = clock64();//end of necessary kernel instructions
//
//		sum1 +=  t2 - t1;
//		if (iter1 == pairs - 1){
//				clockTable[0] = sum1;
//		}
//		if (iter1 == 2 * pairs - 1){
//				clockTable[1] = sum1;
//		}
	}
}


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
	unsigned thread_num = 1;
	long long unsigned data_num = pairs * thread_num;

	////////constant variables
	cuda_mpz_t h_n;
	cuda_mpz_t h_n_;
	cuda_mpz_t h_r2;

	///////get n
	//char n_input[] = "00000038f6e8cfba55dd0e47";
	char n_input[] = "00000003412791fa847ccd00ad83efcae8820aad5457cbd253bc866b3a85184f249ae3a825c6c49af5ebf13cd2ef39ed46a5a0468b153e8521cd5f250049c5491d4f49462edbad1bedb4b48b67f7b59cdb683e6412d40d0000f6e07ba46c0c34d84790e3c83e076c70d3e3eb72ac583700a7664f0efcf67ae4b32254d9d50566357d635b";
	cuda_mpz_set_str_host(&h_n, n_input);

	///////get n_
	//char n__input[] = "0000002e8457440e0d93c489";
	char n__input[] = "00000000f8a46a29539787c065dfe90891abaf64d14a94038141e1a3e0a1712160ad261030b0e23fd01f4d261f6595105cac0cfc6c64730e1992cf8b9403905769f59c60eaa3b2bcd63dc7f616f400be60145b000ec8717ec4fe09d139a2c9d2bf44fbb219fc17f6e6000defe57ce6e46986b33219eb41ef1b4e9df4703a2d658d13eb2d";
	cuda_mpz_set_str_host(&h_n_, n__input);

	///////get r2
	//char r2_input[] = "0000003709d17d8f8686609f";
	char r2_input[] = "000000021433fadceed1a83f846fbf811383c65ef47e61e08788266c019b6b7bcc356c572dc61f969ce683b781633d5f5b3121f23c209d1f6aa8578f32c1c5d8987e6f307c88b649670f861fe820f2288b4164cb3533bd8a70ef57ac229fb6885c20f39eb8acd5f4d79f5b22f9e33b099d08f841cb314a711c390c8d7bcf95b943d69e6b";
	cuda_mpz_set_str_host(&h_r2, r2_input);

	///////get d
	//char d_input[] = "1011011001001001010011110110010101010111001010110101111000111100001";
	char d_input[] = "1011001010001000011110101011010110101110101011010000011101011011100100101110010101101010001111011100010000011011110111011011011101101101100000001000011100011010110010001100110011111000001110111000110010001010001111000001000011110101100011101110011110100100000010000001100001001110101100110111110111010111001000010110100001110110010101111101010110001110010001011111111011101011011111001101010010101001000111111010111011010000011000101101110110000111111011011100011010101010010001101000011001000111110110001101011101101000101011000011010011101001010001011010011100010100101111111001011111111101110001000001001110111111000101110001001001010010100010000101001111111111101110000100000111100001101000010111010000000101011000110110011011110110111111111011001111001110100011101101101011001010100100010101001001010000000110011101110011100000000101100010011101100111101100000011001111010111100000111101111110101110001001000010000111101010110010110100011001111111100001100100000110000110111001101011010000000000000010010010101100000111";
	int d_bitsLength = (int)strlen(d_input);
	int* dBits = (int *) malloc(sizeof(int) * d_bitsLength);

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

	///////get Messages
	long long unsigned mesSize = sizeof(cuda_mpz_t) * data_num;
	cuda_mpz_t *myMes1_h;
	myMes1_h = (cuda_mpz_t*) malloc (mesSize * 2); //CPU, bit1_div and bit0_div lists

	///////time per sample
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc( 2 * sizeof(long long int));	//CPU

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

	///////gmp init
	mpz_t mod;
	mpz_t rand_num;
	mpz_init(mod);
	mpz_init(rand_num);

	mpz_set_str (mod, n_input, 16);

	///////RNG init
	gmp_randstate_t rand_state;
	//gmp_randinit_default (rand_state);
	gmp_randinit_mt(rand_state);
	//gmp_randseed_ui (rand_state, time(NULL));
	gmp_randseed_ui (rand_state, 0);


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

	struct timespec ts1;/////////////////////////////////time
	clock_gettime(CLOCK_REALTIME, &ts1);/////////////////////////////////time

	printf("debug1");

	MontSQMLadder(myMes1_h, pairs, h_r2, h_n, h_n_, dBits, d_bitsLength, clockTable_h);/////////////////////////////////////////kernel

	struct timespec ts2;/////////////////////////////////time
	clock_gettime(CLOCK_REALTIME, &ts2);/////////////////////////////////time
	long long unsigned time_interval = time_diff(ts1, ts2);/////////////////////////////////time
	printf("overall kernel time: %lluns %fms %fs\n", time_interval,  ((double) time_interval) / 1000000,  ((double) time_interval) / 1000000000);/////////////////////////////////time

	///////gmp clear
	gmp_randclear(rand_state);
	mpz_clear(rand_num);
	mpz_clear(mod);

	////////free host
	free(clockTable_h);
	free(myMes1_h);
	free(dBits);

    return 0;
}

