#include <stdio.h>
#include "include/mpz.h"
//
// Nearly minimal CUDA example.
// Compile with:
//
// nvcc -o example example.cu
//

#define N 1000

//
// A function marked __global__
// runs on the GPU but can be called from
// the CPU.
//
// This function multiplies the elements of an array
// of ints by 2.
//
// The entire computation can be thought of as running
// with one thread per array element with blockIdx.x
// identifying the thread.
//
// The comparison i<N is because often it isn't convenient
// to have an exact 1-1 correspondence between threads
// and array elements. Not strictly necessary here.
//
// Note how we're mixing GPU and CPU code in the same source
// file. An alternative way to use CUDA is to keep
// C/C++ code separate from CUDA code and dynamically
// compile and load the CUDA code at runtime, a little
// like how you compile and load OpenGL shaders from
// C/C++ code.
//
__global__
void add(mpz_t* a, mpz_t* b) {
	mpz_t tmp1;
	mpz_init(a);
	mpz_init(b);
	mpz_init(&tmp1);

	mpz_set_i(a,2);
	mpz_set_i(b,3);

	mpz_mult(&tmp1, a, b);
	mpz_set(a, &tmp1);
	//mpz_get_str(&am, addr, N);
    
    //int i = blockIdx.x;
    //if (i<N) {
    //    b[i] = 2*a[i];
    //}
}

__global__
void assign(int* a) {
	*a = 777;
}

int main() {
	// mpz_t a, d;
	// mpz_init(&a);
	// mpz_init(&d);

	// mpz_set_i(a,2);
	// mpz_set_i(b,3);
    //
    // Create int arrays on the CPU.
    // ('h' stands for "host".)
    //
    //int hb[N];

    //
    // Create corresponding int arrays on the GPU.
    // ('d' stands for "device".)
    //
    //char *db;
    //cudaMalloc((void **)&da, N*sizeof(int));
    //cudaMallocManaged((void **)&db, N*sizeof(char));
    mpz_t* a, *b;
    cudaMallocManaged(&a, sizeof(mpz_t));
    cudaMallocManaged(&b, sizeof(mpz_t));

    // int * i;
    // cudaMallocManaged(&i, sizeof(int));
    // printf("i: %d\n", *i);
    // assign<<<1, 1>>>(i);
    // cudaDeviceSynchronize();
    // printf("i: %d\n", *i);


    //
    // Initialise the input data on the CPU.
    //
    // for (int i = 0; i<N; ++i) {
    //     ha[i] = i;
    // }

    //
    // Copy input data to array on GPU.
    //
    //cudaMemcpy(da, ha, N*sizeof(int), cudaMemcpyHostToDevice);

    //
    // Launch GPU code with N threads, one per
    // array element.
    //

    add<<<1, 1>>>(a, b);
    cudaDeviceSynchronize();
    char buf[N];
    //mpz_init(a);
    //mpz_set_i(a,2);
    mpz_get_str(a, buf, N);
    //sprintf(buf, "%d", 100);
    printf("Hello: %s\n", buf);

    

    //
    // Copy output array from GPU back to CPU.
    //
    //cudaMemcpy(hb, db, N*sizeof(char), cudaMemcpyDeviceToHost);

    //char buf [N];

    //mpz_get_str(a, buf, N);
    //printf("%s\n", buf);

    mpz_set_i(a,2);
    mpz_print(a);

    /*for (int i = 0; i<1; ++i) {
        printf("%hhx\n", hb[i]);
    }*/

    //
    // Free up the arrays on the GPU.
    //
    //cudaFree(da);
    cudaFree(&a);
    cudaFree(&b);

    return 0;
}
