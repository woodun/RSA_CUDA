/**
 * @file mpz.c
 *
 * @brief Multiple Precision arithmetic code.
 *
 * @author David Matlack (dmatlack)
 */
#ifndef __418_MPZ_H__
#define __418_MPZ_H__

#include "compile.h"
#include "cuda_string.h"
#include <string.h> //changes

#define LOG2_DIGIT_BASE     32
#define DIGIT_BASE          ((unsigned long long) 1 << (LOG2_DIGIT_BASE))
#define DIGITS_CAPACITY     8 //changes: make enough space for large input
#define MOD_DIGIT_BASE      0xffffffff//changes
#define MOD_LOG2_DIGIT_BASE     31 //changes
#define LOG2_LOG2_DIGIT_BASE 5 //changes
//#define RL 70

#define MPZ_NEGATIVE      1
#define MPZ_NONNEGATIVE  0

typedef unsigned digit_t;

typedef struct {
  digit_t  digits[DIGITS_CAPACITY];
  unsigned words;
  unsigned bits;
} mpz_t;

__device__ __host__ inline void mpz_init(mpz_t *mpz) {
  for (int i = 0; i < DIGITS_CAPACITY; i++) mpz->digits[i] = 0;
  mpz->words = 0;
  mpz->bits = 0;
}

__device__ __host__ inline void mpz_set(mpz_t *to, mpz_t *from) {
  unsigned i;

  #pragma unroll
  for (i = 0; i < from->words; i++) {// changes
    to->digits[i] = from->digits[i];
  }

  #pragma unroll
  for (; i < to->words; i++) {// changes
    to->digits[i] = 0;
  }

  to->words = from->words;
  to->bits = from->bits;
  //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;
}

__host__ inline void mpz_set_str_host(mpz_t *mpz, const char *user_str) {//changes
  unsigned num_digits;
  unsigned i;
  int is_zero;

  #pragma unroll
  for (i = 0; i < DIGITS_CAPACITY; i++) mpz->digits[i] = 0;//changes

  const int bufsize = 1024;
  char buf[bufsize];
  memcpy(buf, user_str, strlen(user_str) + 1);
  buf[bufsize - 1] = (char) 0;
  char *str = &buf[0];

  int len = strlen(str);
  int char_per_digit = LOG2_DIGIT_BASE / 4;
  num_digits = (len + char_per_digit - 1) / char_per_digit;

  is_zero = true;
  #pragma unroll
  for (i = 0; i < num_digits; i ++) {
    str[len - i * char_per_digit] = (char) 0;
    char *start = str + (int) max(len - (i + 1) * char_per_digit, 0);
    digit_t d = strtol(start, NULL, 16);

    /* keep track of whether or not every digit is zero */
    is_zero = is_zero && (d == 0);

    /* parse the string backwards (little endian order) */
    mpz->digits[i] = d;
  }

  mpz->words = num_digits;
  //finding the msb
  digit_t v = mpz->digits[num_digits - 1];
  int msb = 0;

  while (v >>= 1) {
	  msb++;
  }

  mpz->bits = (num_digits - 1) * LOG2_DIGIT_BASE + msb + 1;
  //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;
}

__device__ __host__ inline digit_t digits_add_across(digit_t *digits, unsigned num_digits, digit_t carry) {
  unsigned i = 0;
  unsigned long long value;

  #pragma unroll
  while (carry != 0 && i < num_digits) {
    value = ((unsigned long long) digits[i]) + ((unsigned long long) carry);
    carry  = (digit_t) (value >> LOG2_DIGIT_BASE);
    digits[i] = (digit_t) (value & MOD_DIGIT_BASE);
    i++;
  }

  return carry;
}

__device__ __host__ inline void mpz_mult(mpz_t *dst, mpz_t *op1, mpz_t *op2) {
  unsigned capacity = op1->words + op2->words;

  ///////////////////////debug
  printf("mult:\n");
  printf("dst: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", dst->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", dst->words);
  printf("bits: %u\n", dst->bits);

  printf("op1: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", op1->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", op1->words);
  printf("bits: %u\n", op1->bits);

  printf("op2: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", op2->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", op2->words);
  printf("bits: %u\n", op2->bits);
  printf("##############################################################\n");
  ///////////////////////debug


  #pragma unroll
  for (int i = 0; i < dst->words; i++) {
  //for (int i = word_count; i < DIGITS_CAPACITY; i ++) {
  	  dst->digits[i] = 0;
  }

  digit_t carry;
  digit_t prod;
  unsigned k;
  unsigned long long value;

  digit_t op2_val;
  digit_t op1_val;

  #pragma unroll
  for (unsigned i = 0; i < op2->words; i++) {
	op2_val = op2->digits[i];
	//op2_val = (i < op2->words) ? op2->digits[i] : 0;/////avoiding loads
	#pragma unroll
    for (unsigned j = 0; j < op1->words; j++) {
      k = i + j;
      carry = 0;

      op1_val = op1->digits[j];
      //op1_val = (j < op1->words) ? op1->digits[j] : 0;/////avoiding loads

      value = ((unsigned long long) op2_val) * ((unsigned long long) op1_val) + ((unsigned long long) carry);
      carry  = (digit_t) (value >> LOG2_DIGIT_BASE);
      prod = (digit_t) (value & MOD_DIGIT_BASE);

      digits_add_across(dst->digits + k,     capacity - k,     prod);
      digits_add_across(dst->digits + k + 1, capacity - k - 1, carry);
    }
  }

  unsigned word_count = op1->words + op2->words - 1;
  if(dst->digits[word_count] != 0){
	  word_count++;
  }

  if(dst->digits[word_count - 1] == 0 || dst->digits[word_count] != 0){//check
	  printf("error1!\n");
  }

  unsigned total_bit_count = op1->bits + op2->bits;
  unsigned top_bit_count = ((total_bit_count - 1) & MOD_LOG2_DIGIT_BASE) + 1;
  if( ( dst->digits[word_count - 1] >> (top_bit_count - 1) ) == 0){
	  total_bit_count--;
  }

  top_bit_count = ((total_bit_count - 1) & MOD_LOG2_DIGIT_BASE) + 1;
  if( ( dst->digits[word_count - 1] >> (top_bit_count - 1) ) != 1 ){
	  printf("error2! %x\n", ( dst->digits[word_count - 1] >> (top_bit_count - 1) ) );//check
  }

  dst->words = word_count;
  dst->bits = total_bit_count;


  ///////////////////////debug
  printf("mult:\n");
  printf("dst: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", dst->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", dst->words);
  printf("bits: %u\n", dst->bits);

  printf("op1: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", op1->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", op1->words);
  printf("bits: %u\n", op1->bits);

  printf("op2: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", op2->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", op2->words);
  printf("bits: %u\n", op2->bits);
  printf("##############################################################\n");
  printf("##############################################################\n");
  ///////////////////////debug
  //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;
}

__device__ __host__ inline void mpz_bitwise_truncate(mpz_t *dst, mpz_t *src, int RL) {//changes

    ///////////////////////debug
    printf("truncate:\n");
    printf("dst: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", dst->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", dst->words);
    printf("bits: %u\n", dst->bits);

    printf("src: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", src->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", src->words);
    printf("bits: %u\n", src->bits);
    printf("##############################################################\n");
    ///////////////////////debug

//  if(RL >= dst->bits){
//	  return;
//  }

  digit_t *src_digits = src->digits;
  digit_t *dst_digits = dst->digits;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  //int ls_digits = capacity - rs_digits - 1;
  int ls_remainder = LOG2_DIGIT_BASE - ( RL & MOD_LOG2_DIGIT_BASE );
  //int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  #pragma unroll
  for(int d_index = dst->words - 1; d_index > rs_digits; d_index--) {
	  dst_digits[d_index] = 0;
  }

  digit_t top_word = src_digits[rs_digits] & ( 0xffffffff >> ls_remainder);
  dst_digits[rs_digits] = top_word;

  #pragma unroll
  for(int d_index = rs_digits - 1; d_index >= 0; d_index--) {
	  dst_digits[d_index] = src_digits[d_index];
  }

  unsigned word_count = rs_digits + 1;
  unsigned total_bit_count = RL;
  unsigned top_bit_count = ((total_bit_count - 1) & MOD_LOG2_DIGIT_BASE) + 1;/////////////////////check when rl = 64

  //finding the msb, for efficiency we assume heading zeros does not pass word boundary
   while ( (top_word >> ( top_bit_count - 1 ) ) == 0 ) {
	   top_bit_count--;
  }

  if( (top_word >> ( top_bit_count - 1 ) ) != 1 ){
	  printf("error6! %x\n", (top_word >> ( top_bit_count - 1 ) ) );//check
  }

  dst->bits = (word_count - 1) * LOG2_DIGIT_BASE + top_bit_count;
  dst->words = word_count;

  ///////////////////////debug
  printf("truncate:\n");
  printf("dst: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", dst->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", dst->words);
  printf("bits: %u\n", dst->bits);

  printf("src: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", src->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", src->words);
  printf("bits: %u\n", src->bits);
  printf("##############################################################\n");
  printf("##############################################################\n");
  ///////////////////////debug
}

__device__ __host__ inline void mpz_bitwise_truncate_eq(mpz_t *mpz, int RL) {//changes

    ///////////////////////debug
    printf("truncateeq:\n");
    printf("mpz: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", mpz->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", mpz->words);
    printf("bits: %u\n", mpz->bits);
    printf("##############################################################\n");
    ///////////////////////debug

//  if(RL >= dst->bits){
//	  return;
//  }

  digit_t *digits = mpz->digits;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  //int ls_digits = capacity - rs_digits - 1;
  int ls_remainder = LOG2_DIGIT_BASE - ( RL & MOD_LOG2_DIGIT_BASE );
  //int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  #pragma unroll
  for(int d_index = mpz->words - 1; d_index > rs_digits; d_index--) {//constant time for specific rl
	digits[d_index] = 0;
  }

  digit_t top_word = digits[rs_digits] & ( 0xffffffff >> ls_remainder);
  digits[rs_digits] = top_word;

  mpz->bits = RL;
  mpz->words = (RL + LOG2_DIGIT_BASE - 1 ) >> LOG2_LOG2_DIGIT_BASE;

  unsigned word_count = rs_digits + 1;
  unsigned total_bit_count = RL;
  unsigned top_bit_count = ((total_bit_count - 1) & MOD_LOG2_DIGIT_BASE) + 1;/////////////////////check when rl = 64

  //finding the msb, for efficiency we assume heading zeros does not pass word boundary
   while ( (top_word >> ( top_bit_count - 1 ) ) == 0 ) {
	   top_bit_count--;
  }

  if( (top_word >> ( top_bit_count - 1 ) ) != 1 ){
	  printf("error6! %x\n", (top_word >> ( top_bit_count - 1 ) ) );//check
  }

  mpz->bits = (word_count - 1) * LOG2_DIGIT_BASE + top_bit_count;
  mpz->words = word_count;

  ///////////////////////debug
  printf("truncateeq:\n");
  printf("mpz: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", mpz->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", mpz->words);
  printf("bits: %u\n", mpz->bits);
  printf("##############################################################\n");
  printf("##############################################################\n");
  ///////////////////////debug
}

__device__ __host__ inline int mpz_compare(mpz_t *a, mpz_t *b) {

  if(a->bits > b->bits){
	  return 1;
  }
  if(a->bits < b->bits){
	  return -1;
  }
  unsigned top_index = a->words - 1;
  if(a->digits[top_index] > b->digits[top_index]){
	  return 1;
  }
  if(a->digits[top_index] < b->digits[top_index]){
	  return -1;
  }
  return 0;
}

__device__ __host__ inline int mpz_gte(mpz_t *a, mpz_t *b) {
  return (mpz_compare(a, b) >= 0);
}

__device__ __host__ inline void mpz_bitwise_rshift_eq(mpz_t *mpz, int RL) {//changes

//  if(RL >= mpz->bits){
//	  for (int i = 0; i < mpz->words; i++) mpz->digits[i] = 0;
//	  mpz->words = 0;
//	  mpz->bits = 0;
//	  return;
//  }

  digit_t *digits = mpz->digits;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  int rs_remainder = RL & MOD_LOG2_DIGIT_BASE;

  #pragma unroll
  for(int d_index = 0; d_index < mpz->words - 1 - rs_digits; d_index++) {
	digits[d_index] = ( digits[d_index + rs_digits] >> rs_remainder ) | ( digits[d_index + rs_digits + 1] << ( LOG2_DIGIT_BASE - rs_remainder ) );
  }

  digits[mpz->words - 1 - rs_digits] = digits[mpz->words - 1] >> rs_remainder;

  #pragma unroll
  for(int d_index = mpz->words - rs_digits; d_index <= mpz->words - 1; d_index++) {//could cause overflow, but not in RSA
  	digits[d_index] = 0;
  }

  mpz->bits = mpz->bits - RL;
  mpz->words = (mpz->bits + LOG2_DIGIT_BASE - 1 ) >> LOG2_LOG2_DIGIT_BASE;
}

__device__ __host__ inline void mpz_bitwise_rshift(mpz_t *dst, mpz_t *src, int RL) {//changes

    ///////////////////////debug
    printf("rshift:\n");
    printf("dst: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", dst->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", dst->words);
    printf("bits: %u\n", dst->bits);

    printf("src: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", src->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", src->words);
    printf("bits: %u\n", src->bits);
    printf("##############################################################\n");
    ///////////////////////debug

//  if(RL >= dst->bits){
//	  for (int i = 0; i < dst->words; i++) dst->digits[i] = 0;
//	  dst->words = 0;
//	  dst->bits = 0;
//	  return;
//  }

  digit_t *src_digits = src->digits;
  digit_t *dst_digits = dst->digits;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  int rs_remainder = RL & MOD_LOG2_DIGIT_BASE;

  #pragma unroll
  for(int d_index = 0; d_index < src->words - 1 - rs_digits; d_index++) {
	  dst_digits[d_index] = ( src_digits[d_index + rs_digits] >> rs_remainder ) | ( src_digits[d_index + rs_digits + 1] << ( LOG2_DIGIT_BASE - rs_remainder ) );
  }

  dst_digits[src->words - 1 - rs_digits] = src_digits[src->words - 1] >> rs_remainder;

  #pragma unroll
  for(int d_index = src->words - rs_digits; d_index <= dst->words - 1; d_index++) {//could cause overflow, but not in RSA
	  dst_digits[d_index] = 0;
  }

  dst->bits = src->bits - RL;
  dst->words = (dst->bits + LOG2_DIGIT_BASE - 1 ) >> LOG2_LOG2_DIGIT_BASE;


  ///////////////////////debug
  printf("rshift:\n");
  printf("dst: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", dst->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", dst->words);
  printf("bits: %u\n", dst->bits);

  printf("src: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", src->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", src->words);
  printf("bits: %u\n", src->bits);
  printf("##############################################################\n");
  printf("##############################################################\n");
  ///////////////////////debug
}

__device__ __host__ inline void mpz_add(mpz_t *dst, mpz_t *op1, mpz_t *op2) {

	  ///////////////////////debug
	  printf("add:\n");
	  printf("dst: \n");
	  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
		  printf("%08x", dst->digits[i]);
	  }
	  printf("\n");
	  printf("words: %u\n", dst->words);
	  printf("bits: %u\n", dst->bits);

	  printf("op1: \n");
	  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
		  printf("%08x", op1->digits[i]);
	  }
	  printf("\n");
	  printf("words: %u\n", op1->words);
	  printf("bits: %u\n", op1->bits);

	  printf("op2: \n");
	  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
		  printf("%08x", op2->digits[i]);
	  }
	  printf("\n");
	  printf("words: %u\n", op2->words);
	  printf("bits: %u\n", op2->bits);
	  printf("##############################################################\n");
	  ///////////////////////debug

  unsigned capacity = max(op1->words, op2->words);

  digit_t carry = 0;
  digit_t a;
  digit_t b;
  unsigned long long value;

  for (unsigned i = 0; i < capacity + 1; i++) {
    a = (i < op1->words) ? op1->digits[i] : 0;
    b = (i < op2->words) ? op2->digits[i] : 0;

    value = ((unsigned long long) a) + ((unsigned long long) b) + ((unsigned long long) carry);
    carry  = (digit_t) (value >> LOG2_DIGIT_BASE);
    dst->digits[i] = (digit_t) (value & MOD_DIGIT_BASE);
  }

  unsigned word_count;
  if(value != 0){
	  word_count = capacity + 1;
  }else{
	  word_count = capacity;
  }

  #pragma unroll
  for (int i = word_count; i < dst->words; i++) {
  //for (int i = word_count; i < DIGITS_CAPACITY; i ++) {
  	  dst->digits[i] = 0;
  }

  unsigned total_bit_count = max(op1->bits, op2->bits);
  unsigned top_bit_count = ((total_bit_count - 1) & MOD_LOG2_DIGIT_BASE) + 1;
  if( ( dst->digits[word_count - 1] >> (top_bit_count - 1) ) != 1){
	  total_bit_count++;
  }

  top_bit_count = ((total_bit_count - 1) & MOD_LOG2_DIGIT_BASE) + 1;
  if(dst->digits[word_count - 1] == 0 || dst->digits[word_count] != 0){//check
	  printf("error3!\n");
  }
  if( ( dst->digits[word_count - 1] >> (top_bit_count - 1) ) != 1 ){
	  printf("error4! %x\n", ( dst->digits[word_count - 1] >> (top_bit_count - 1) ) );//check
  }

  dst->words = word_count;
  dst->bits = total_bit_count;
  //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;


  ///////////////////////debug
  printf("add:\n");
  printf("dst: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", dst->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", dst->words);
  printf("bits: %u\n", dst->bits);

  printf("op1: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", op1->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", op1->words);
  printf("bits: %u\n", op1->bits);

  printf("op2: \n");
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
	  printf("%08x", op2->digits[i]);
  }
  printf("\n");
  printf("words: %u\n", op2->words);
  printf("bits: %u\n", op2->bits);
  printf("##############################################################\n");
  printf("##############################################################\n");
  ///////////////////////debug
}

__device__ __host__ inline void mpz_sub(mpz_t *dst, mpz_t *op1, mpz_t *op2) {

	  ///////////////////////debug
	  printf("sub:\n");
	  printf("dst: \n");
	  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
		  printf("%08x", dst->digits[i]);
	  }
	  printf("\n");
	  printf("words: %u\n", dst->words);
	  printf("bits: %u\n", dst->bits);

	  printf("op1: \n");
	  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
		  printf("%08x", op1->digits[i]);
	  }
	  printf("\n");
	  printf("words: %u\n", op1->words);
	  printf("bits: %u\n", op1->bits);

	  printf("op2: \n");
	  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
		  printf("%08x", op2->digits[i]);
	  }
	  printf("\n");
	  printf("words: %u\n", op2->words);
	  printf("bits: %u\n", op2->bits);
	  printf("##############################################################\n");
	  ///////////////////////debug


	unsigned capacity = max(op1->words, op2->words);

	// Complement each digit by subtracting it from BASE-1
	#pragma unroll
	for (unsigned i = 0; i < capacity; i++) {
		dst->digits[i] = (digit_t) ((DIGIT_BASE - 1) - op2->digits[i]);
	}

	// Add 1
	digits_add_across(dst->digits, capacity, 1);

	#pragma unroll
	for (int i = capacity; i < dst->words; i++) {
	//for (int i = word_count; i < DIGITS_CAPACITY; i ++) {
		  dst->digits[i] = 0;
	}

    digit_t carry = 0;
    digit_t a;
    digit_t b;
    unsigned long long value;
    digit_t tmp;
    digit_t top_digit = 0;
    int word_count = -1;

    for (unsigned i = 0; i < capacity; i++) {
      //a = (i < op1->words) ? op1->digits[i] : 0;
      //b = (i < op2->words) ? dst->digits[i] : 0;
       a = op1->digits[i];
       b = dst->digits[i];

      value = ((unsigned long long) a) + ((unsigned long long) b) + ((unsigned long long) carry);
      carry  = (digit_t) (value >> LOG2_DIGIT_BASE);

      tmp = (digit_t) (value & MOD_DIGIT_BASE);
      dst->digits[i] = tmp;
      if(tmp != 0){
    	  word_count = i;
    	  top_digit = tmp;
      }
    }

    if(carry != 1){////////in RSA it must be positive number
    	printf("error5! negative!\n");//check
    }

    word_count++;
    dst->words = word_count;

//    if(word_count == 0){//////should not happen in RSA
//    	dst->bits = 0;
//    	return;
//    }

    //finding the msb
    int msb = 0;
    while (top_digit >>= 1) {
  	  msb++;
    }

    dst->bits = (word_count - 1) * LOG2_DIGIT_BASE + msb + 1;
    //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;

    ///////////////////////debug
    printf("sub:\n");
    printf("dst: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", dst->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", dst->words);
    printf("bits: %u\n", dst->bits);

    printf("op1: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", op1->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", op1->words);
    printf("bits: %u\n", op1->bits);

    printf("op2: \n");
    for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
  	  printf("%08x", op2->digits[i]);
    }
    printf("\n");
    printf("words: %u\n", op2->words);
    printf("bits: %u\n", op2->bits);
    printf("##############################################################\n");
    printf("##############################################################\n");
    ///////////////////////debug
}

__device__ __host__ inline digit_t mpz_get_last_digit(mpz_t *mpz) {//changes
	return mpz->digits[0];
}

__host__ inline char* mpz_get_str(mpz_t *mpz, char *str, int bufsize) {
  int print_zeroes = 0; // don't print leading 0s
  int i;
  int str_index = 0;

  #pragma unroll
  for (i = DIGITS_CAPACITY - 1; i >= 0; i--) {
    unsigned digit = mpz->digits[i];

    if (digit != 0 || print_zeroes) {
      if (bufsize < str_index + 8) {
        return NULL;
      }
      if (!print_zeroes) {
        str_index += sprintf(str + str_index, "%x", digit);
      }
      else {
        str_index += sprintf(str + str_index, "%08x", digit);
      }
      print_zeroes = 1;
    }
  }

  str[str_index] = (char) 0;

  /* the number is zero */
  if (print_zeroes == 0) {
    str[0] = '0';
    str[1] = (char) 0;
  }

  return str;
}

__device__ inline void mpz_print_str_device(mpz_t *mpz) {//changes
  int print_zeroes = 0; // don't print leading 0s

  #pragma unroll
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
    unsigned digit = mpz->digits[i];

    if (digit != 0 || print_zeroes) {
      if (!print_zeroes) {
    	  printf("%x", digit);
      }
      else {
    	  printf("%08x", digit);
      }
      print_zeroes = 1;
    }
  }
  /* the number is zero */
  if (print_zeroes == 0) {
	  printf("0");
  }
}


#endif /* __418_MPZ_H__ */
