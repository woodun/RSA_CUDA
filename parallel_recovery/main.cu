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
	if (argc < 3){
		printf("device id required.\n");
		exit(EXIT_FAILURE);
	}

	int peak_clk = 1;//kHz
	int dev_id = 0;
	cudaSetDevice(dev_id);
	cudaDeviceGetAttribute(&peak_clk, cudaDevAttrClockRate, dev_id);
	float clock_rate = (float) peak_clk;
	printf("clock_rate_out_kernel:%f\n", clock_rate);

	long devid = strtol(argv[2], NULL, 10);
	cudaSetDevice(devid);

	long x = strtol(argv[1], NULL, 10);
	long long unsigned pairs = x;
	unsigned thread_num = 2;
	long long unsigned data_num = pairs * thread_num;

	////////constant variables
	cuda_mpz_t h_n;
	cuda_mpz_t h_n_;
	cuda_mpz_t h_r2;

	///////get n
	//char n_input[] = "00000038f6e8cfba55dd0e47";
	//char n_input[] = "00000003412791fa847ccd00ad83efcae8820aad5457cbd253bc866b3a85184f249ae3a825c6c49af5ebf13cd2ef39ed46a5a0468b153e8521cd5f250049c5491d4f49462edbad1bedb4b48b67f7b59cdb683e6412d40d0000f6e07ba46c0c34d84790e3c83e076c70d3e3eb72ac583700a7664f0efcf67ae4b32254d9d50566357d635b";
	char n_input[] = "000000017eafe8bfb746bbf8a18e6febb4b3bd2fe5ff65adc76b1b65472008e722ae2d436d03997f4e6e16d101a6bb15c2f85c285b8fd03c92f5d62975e790dec7b58e58c012c0dbd9dfa26dd8e9191b54f7f9e8503678d60fd97555334b35eaeba92783edc8f3a976c2c5a7f288f14f819622b8989f56873f4fcb0e6ffdef89d6c52dd93b312aaaaca520e396b01bbb8983a14e42fd5f2c3551af68bab481614857970f8646a6558c45212185da9bfe9ee743e91f7f99439b4ecc1e4756a849950a200fb31a0f3e46c7b642af57194a79bd1aa486b08862ac199fee1b96ec2f5bdac157c935a2b405693b78bb3b55afbcf288750bf43b539e5dadf2a2a1d08bdfc6d577";
	cuda_mpz_set_str_host(&h_n, n_input);

	///////get n_
	//char n__input[] = "0000002e8457440e0d93c489";
	//char n__input[] = "00000000f8a46a29539787c065dfe90891abaf64d14a94038141e1a3e0a1712160ad261030b0e23fd01f4d261f6595105cac0cfc6c64730e1992cf8b9403905769f59c60eaa3b2bcd63dc7f616f400be60145b000ec8717ec4fe09d139a2c9d2bf44fbb219fc17f6e6000defe57ce6e46986b33219eb41ef1b4e9df4703a2d658d13eb2d";
	char n__input[] = "000000001edc182fc42caa0bd519f16df49b672f78b08d904aeba23553a9b4739e3c8ac75abbc7c1f4eb3fd3d008b97fcba580667737e9f5c5a54a3ed0d4bdcf802c04c035575dcd1d4ddf980f368a7fd4cd8e0d574cfb97ebce7fec7a199c3af5928918b69c2b0ace65b9715c11ea379a8bc3ee72c2507158d64f6991bea04bac03da924bdb7fbcfd173133cf066e3abe62324689c2b7b59af7fc3e687737017b0e3b4db2c35e52378d3af356b5fff18bafc287d9ae5a3ad028f5358993c2e62bc87c4c89d2059a523951b6ca196ee35bb813bc6c1f78b17843367eeb0bce306790284b2ad8ff3b4a47e69451cc27b180f75b5dd2253269eadb87a9f9f6cd5914616bb9";
	cuda_mpz_set_str_host(&h_n_, n__input);

	///////get r2
	//char r2_input[] = "0000003709d17d8f8686609f";
	//char r2_input[] = "000000021433fadceed1a83f846fbf811383c65ef47e61e08788266c019b6b7bcc356c572dc61f969ce683b781633d5f5b3121f23c209d1f6aa8578f32c1c5d8987e6f307c88b649670f861fe820f2288b4164cb3533bd8a70ef57ac229fb6885c20f39eb8acd5f4d79f5b22f9e33b099d08f841cb314a711c390c8d7bcf95b943d69e6b";
	char r2_input[] = "000000008abb423bfadf1233fc46527546f44c092476681a687cc1845a29e71d5f2a0a912b83dd51a5cdd8f1d577ebef35e5bc7b40aeae83abdd4c892e01222e5e610a318b6b4457bd145b505c01d153de75ebb1996d2001a55d9ccabc2c61d7aac2509d60033bc3302c783841df3f3f3a264a1284d57cdd4b2b7f91cc3b783a89264bc435c285f77e66035476f100602cd48b6d7ff4d1acb21704096561cc73b1c71caebbe6c14f073a42a99d67eb0667d73e5d88f1f860dc9796ada3eb78c2373ce1321bb506c07789a1fdbad607a589f547e54e97e048e1a35aad06d75cef88ef2575c77d4550d8894c8c135372663b2d783c682aae6ea9bc66d8f326542293665a99";
	cuda_mpz_set_str_host(&h_r2, r2_input);

	///////get d
	//char d_input[] = "1011011001001001010011110110010101010111001010110101111000111100001";
	//char d_input[] = "1011001010001000011110101011010110101110101011010000011101011011100100101110010101101010001111011100010000011011110111011011011101101101100000001000011100011010110010001100110011111000001110111000110010001010001111000001000011110101100011101110011110100100000010000001100001001110101100110111110111010111001000010110100001110110010101111101010110001110010001011111111011101011011111001101010010101001000111111010111011010000011000101101110110000111111011011100011010101010010001101000011001000111110110001101011101101000101011000011010011101001010001011010011100010100101111111001011111111101110001000001001110111111000101110001001001010010100010000101001111111111101110000100000111100001101000010111010000000101011000110110011011110110111111111011001111001110100011101101101011001010100100010101001001010000000110011101110011100000000101100010011101100111101100000011001111010111100000111101111110101110001001000010000111101010110010110100011001111111100001100100000110000110111001101011010000000000000010010010101100000111";
	char d_input[] = "10001000101011001000100111111011010100111011110111010101011010110001010100100000100101011010111110011011111101110000110010110101101011011001001000010010000001110011010011101111011001010011011001110100110101001001010101110111000111101010101111101011100110000001010010100101110110110110010001010010110111100010110010111000010111000000010010110000100001111100010110100001110101111100010101000101010001011010010111001100011111011010000011110001000011101100111010101110001000010111010000100010110000001101011101101000110101101110001000100000010011101000010010101011010011000100101111001101011101111101001000011100000011000000111101101011100010011101001110000001001010110010011111100001000101101100111001111001110110110111011001001010000111010000101100000101100011100001110011010100111011000101011100000101101010100110101000100010000001010001111111000011001100011001110001100101001000110101010110001011000100011110111111000011011110010111001000001010001101100011110000000011011011001111101000011110111100010100011001101011110011001011010101000111001010000110001000000110000110010010110111100100101100001000111010100101101100100110110110011111110100001000101110101010111010110010100001100111001010111101011000001000110100100001100111111010010111111000101000010111110010001011000110101010000010110111111101001100000010110010001110000001101011111110100000011101110010010010101010000010111111100001001010000010010111010011100111111000011001001111000111000001010110101011100111000110111110100100010010101101111100011101001101110111111001010001011101000001000100100101010000001000001010000111001111101101111010111101000001111011011011100010010101101100011001100010010101001000000000101101110010110000001001011010100001100111010001100011100101100001011110011111110011100011111010111111011000100101000110101000010100011101101011010010100100100011001001100001001001011000110100110100101101001000010110111101011001010100001100110011110001101001010001111001001101011001110101001010000011000001111111000010110001000000110110010000000011100101011001100111011001010011";
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

	///////get Messages
	long long unsigned mesSize = sizeof(cuda_mpz_t) * data_num;
	cuda_mpz_t *myMes1_h;
	myMes1_h = (cuda_mpz_t*) malloc (mesSize * 2); //CPU, bit1_div and bit0_div lists
	cuda_mpz_t *myMes1_d;
	cudaMalloc((cuda_mpz_t **) &myMes1_d, mesSize * 2); //GPU

	for (long long unsigned p = 0; p < 2 * data_num; ++p){
		cuda_mpz_init( &myMes1_h[p]);
	}

	///////time per sample
	long long int *clockTable_h;
	clockTable_h = (long long int*) malloc( 2 * pairs * sizeof(long long int));	//CPU
	long long int *clockTable_d;
	cudaMalloc((void **) &clockTable_d, 2 * pairs * sizeof(long long int)); //GPU

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
//	known_bits[1] = 0;
//	known_bits[2] = 1;
	int known_bits_length = 1;
	int div_con = 0;
	int wrong_key = 0;

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
	gmp_randseed_ui (rand_state, time(NULL));
	//gmp_randseed_ui (rand_state, 0);

	printf("current bits: ");
	for(int i = 0; i < known_bits_length; i++){
		printf("%d", known_bits[i]);
	}
	printf("\n");

	int vote1 = 0;
	int vote0 = 0;
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

			if (div_con == 1 && bit1_div_num < 2 * data_num){//0,1,4,5,...
				cuda_mpz_set( &myMes1_h[bit1_div_num], &r1);
				bit1_div_num++;
				cuda_mpz_set( &myMes1_h[bit1_div_num], &r2);
				bit1_div_num++;
				bit1_div_num+=2;
			}
			if (div_con == 4 && bit0_div_num < 2 * data_num){//2,3,6,7,...
				bit0_div_num+=2;
				cuda_mpz_set( &myMes1_h[bit0_div_num], &r1);
				bit0_div_num++;
				cuda_mpz_set( &myMes1_h[bit0_div_num], &r2);
				bit0_div_num++;
			}
			if (bit1_div_num == 2 * data_num && bit0_div_num == 2 * data_num){
				break;
			}
		}

		cudaMemcpy(myMes1_d, myMes1_h, mesSize * 2 , cudaMemcpyHostToDevice);///////////////bit1_div and bit0_div lists

		cudaFuncSetAttribute(MontSQMLadder, cudaFuncAttributePreferredSharedMemoryCarveout, 100);

		struct timespec ts1;/////////////////////////////////time
		clock_gettime(CLOCK_REALTIME, &ts1);/////////////////////////////////time

		MontSQMLadder<<<2 * pairs, 2>>>(myMes1_d, h_r2, h_n, h_n_, dBits_d, d_bitsLength, clockTable_d);/////////////////////////////////////////kernel
		cudaDeviceSynchronize();

		struct timespec ts2;/////////////////////////////////time
		clock_gettime(CLOCK_REALTIME, &ts2);/////////////////////////////////time
		long long unsigned time_interval = time_diff(ts1, ts2);/////////////////////////////////time
		printf("overall kernel time: %lluns %fms %fs\n", time_interval,  ((double) time_interval) / 1000000,  ((double) time_interval) / 1000000000);/////////////////////////////////time

		cudaMemcpy(clockTable_h, clockTable_d, 2 * pairs * sizeof(long long int), cudaMemcpyDeviceToHost);

		double sum_time1 = 0;
		for (long long unsigned q = 0; q < 2 * pairs; q+=2){
			//fprintf(fp1, "%d ", clockTable_h[q]);///////div count///////////////////print file
			sum_time1 += clockTable_h[q];
		}
		sum_time1 = sum_time1 / pairs;

		double sum_time4 = 0;
		for (long long unsigned q = 1; q < 2 * pairs; q+=2){
			//fprintf(fp1, "%d ", clockTable_h[q]);///////div count///////////////////print file
			sum_time4 += clockTable_h[q];
		}
		sum_time4 = sum_time4 / pairs;

		double diff = sum_time1 - sum_time4;
	//	printf("%f %f %f\n", sum_time1, sum_time4, diff);
	//	printf("%f\n", diff);
		printf ("bit1_div: %fms %fcycles ", sum_time1 / clock_rate, sum_time1);
		printf ("bit0_div: %fms %fcycles ", sum_time4 / clock_rate, sum_time4);
		printf ("difference: %fms %fcycles\n", diff / clock_rate, diff);

		if(diff > 30000){//bit is 1
			known_bits[known_bits_length] = 1;
			printf("bit is 1.\n");
		}else if(diff < -30000){//bit is 0
			known_bits[known_bits_length] = 0;
			printf("bit is 0.\n");
		}else{//EOB
			//printf("end of bits.\n");

			if(diff > 10000){//bit is 1
				vote1++;
				printf("vote 1.\n");
			}else if(diff < -10000){//bit is 0
				vote0++;
				printf("vote 0.\n");
			}else{
				printf("result is discarded.\n");
				continue;
			}

			if( vote1 >= 3 ){/////////////////////////////////if not accepted for too many times, then decide by voting
				known_bits[known_bits_length] = 1;
				printf("bit is voted 1.\n");
			}else if( vote0 >= 3 ){
				known_bits[known_bits_length] = 0;
				printf("bit is voted 0.\n");
			}else{
				printf("bit not acceptable.\n");
				continue;
			}
		}
		vote1 = 0;
		vote0 = 0;

		known_bits_length++;

		printf("current bits: ");
		for(int i = 0; i < known_bits_length; i++){
			printf("%d", known_bits[i]);
		}
		printf("\n");

		if(known_bits[known_bits_length - 1] != dBits[known_bits_length - 1]){
			wrong_key = 1;
			printf("wrong key!\n");
			break;
		}
		fflush(stdout);
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

	printf("pair count: %llu\n", pairs);
	fflush(stdout);

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

