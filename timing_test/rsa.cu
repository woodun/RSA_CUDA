//#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include "include/mpz.h"
#include <assert.h>
#include <unistd.h>
//static char usage[] = "usage: %s -n n_samles -t n_div_thread_limit\n";
#define N 32

// find R
// __device__ __host__ int find_r(mpz * r, mpz * n){
//     for(int i = DIGITS_CAPACITY * LOG2_DIGIT_BASE; i>=0; i--){

//     }
// }

// __global__
// void rsa(uint64_t* ms, int n, long long int* t) {
//     long long int t1, t2;
//     uint64_t i;
//     //uint64_t temp;
//     t[threadIdx.x] = threadIdx.x;
//     a[threadIdx.x] = 0;
//     for(i=0;i<10000000;i++){
//         a[threadIdx.x] += 88;
//         if(n==3333333){ //never true, to make sure line above is executed
//             a[threadIdx.x] = 11111;//powf(float(a[threadIdx.x]),float(n));
//         }
//     }
//     t1 = clock64();
//     //temp = powf(float(a[threadIdx.x]),float(n));
//     if(threadIdx.x >= n){ // branch that causes divergence
//         // for(uint64_t j=0;j<10000;j++){
//             //a[threadIdx.x] += powf(float(a[threadIdx.x]),float(n));
            
//         //asm("sub.s64 %rd20, %rd19, 77777;")
//         t[threadIdx.x] += 22;

//         // }
//         //t[threadIdx.x] = 99;
//         //for(j=0;j<77777;j++){
//             // if(n==5555555){ //never true, to make sure line above is executed
//             //     a[threadIdx.x] = 6666;//powf(float(a[threadIdx.x]),float(n));
//             // }
//         //}
//         //a[threadIdx.x] -= temp - 777;
//     } else {
//         //t[threadIdx.x] *= 33;
//         //asm("add.s64 %rd32, %rd33, %rd34;");
//         //a[threadIdx.x] += powf(float(a[threadIdx.x]),float(n));
//         //asm("add.s64 %rd20, %rd19, 77777;");
//     }
//     t2 = clock64();
//     a[threadIdx.x] = t2 - t1; 
//     t[threadIdx.x] = threadIdx.x;
// }



int main(int argc, char ** argv) {
    // char c;
    // //int n_samples = -1;
    // int div_threads = -1;
    // while ((c = getopt(argc, argv, "t:n:")) != -1)
    //   switch (c) {
    //   case 'n':
    //     n_samples = atoi(optarg);
    //     break;
    //   case 't':
    //     div_threads = atoi(optarg);
    //     break;
    //   case '?':
    //     fprintf(stderr, usage, argv[0]);
    //     exit(1);
    //     break;
    // }
    mpz_t *a, *b, *t;
    a = (mpz_t*) malloc(sizeof(mpz_t));
    b = (mpz_t*) malloc(sizeof(mpz_t));
    t = (mpz_t*) malloc(sizeof(mpz_t));
    mpz_init(a);
    mpz_init(b);
    mpz_init(t);
    mpz_set_str(a, "12345679");
    mpz_print(a);
    printf("a is zero? (0): %d\n", mpz_is_zero(a));
    mpz_bit_rshift(a,2);
    mpz_print(a);
    mpz_bit_rshift(a,26);
    mpz_print(a);
    mpz_bit_rshift(a,1000000);
    mpz_print(a);
    printf("a is zero? (1): %d\n", mpz_is_zero(a));

    //printf("a is zero? (0): %d\n", mpz_is_zero(a));
    //mpz_set_str(a, "0");
    //printf("a is zero? (1): %d\n", mpz_is_zero(a));
    //printf("a is zero? (0): %d\n", mpz_is_zero(a));
    //mpz_print(a);
    //mpz_bit_rshift(a);
    //printf("a is zero? (1): %d\n", mpz_is_zero(a));


    printf("The end\n");

    // if(n_samples <= 0 || div_threads <= 0){
    //     fprintf(stderr, "Wrong parameters! Must specify -n and -t\n");
    //     exit(1);
    // }

    //allocate
    // uint64_t *probes;
    // probes = (uint64_t*) malloc(n_samples * sizeof(uint64_t));
    // assert(probes);
    // long long int* a, *t;
    // long long int* a_d, *t_d;
    // cudaMalloc(&a_d, sizeof(long long int) * N);
    // cudaMalloc(&t_d, sizeof(long long int) * N);
    // a = (long long int *) malloc(sizeof(long long int) * N);
    // t = (long long int *) malloc(sizeof(long long int) * N);
    //cudaMemset(a_d, 333, N);
    //cudaEvent_t start, stop;
    //float time;
    //cudaEventCreate(&start);
    //cudaEventCreate(&stop);
    
    // for(int i=0;i<n_samples;i++){
    //     //*a = 1;
    //     //cudaEventRecord(start, 0);
    //     loop<<<1, N>>>(a_d, 31, t_d);
        
    //     //cudaEventRecord(stop, 0);
    //     //cudaEventSynchronize(stop);
    //     cudaDeviceSynchronize();
    //     cudaMemcpy(a, a_d, sizeof(long long int) * N, cudaMemcpyDefault);
    //     cudaMemcpy(t, t_d, sizeof(long long int) * N, cudaMemcpyDefault);

    //     //cudaEventElapsedTime(&time, start, stop);
    //     printf ("%lld\n", *a);

    // }

    // printf("---------\n");

    // for(int i=0;i<n_samples;i++){
    //     //*a = 1;
    //     //cudaEventRecord(start, 0);
    //     loop<<<1, N>>>(a_d, 100, t_d);
        
    //     //cudaEventRecord(stop, 0);
    //     //cudaEventSynchronize(stop);
    //     cudaDeviceSynchronize();
    //     cudaMemcpy(a, a_d, sizeof(long long int) * N, cudaMemcpyDefault);
    //     cudaMemcpy(t, t_d, sizeof(long long int) * N, cudaMemcpyDefault);

    //     //cudaEventElapsedTime(&time, start, stop);
    //     printf ("%lld\n", *a);

    // }

    // // // printf("a: %d\n", *a);
    // // // for (int j=0;j<N;j++){
    // // //     printf("%d:%d\t", j, t[j]);
    // // // }
    // // // printf("\n");


    // cudaFree(&a);
    // cudaFree(&t);
    // free(a);
    // free(t);

    // Retrieve result from device and store it in host array

    return 0;
}
