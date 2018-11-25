/**
 * @file cuda_mpz.c
 *
 * @brief Multiple Precision arithmetic code.
 *
 * @author David Matlack (dmatlack)
 */
#ifndef __418_MPZ_H__
#define __418_MPZ_H__

#include "compile.h"
#include "cuda_string.h"
#include "digit.h"
#include <string.h> //changes
#include <gmp.h>

#define MPZ_NEGATIVE      1
#define MPZ_NONNEGATIVE  0

/** @brief struct used to represent multiple precision integers (Z). */
typedef struct {
  digit_t  digits[DIGITS_CAPACITY];
  unsigned capacity;
  unsigned words;
  unsigned bits;
  char     sign;
} cuda_mpz_t;

__host__ inline char* cuda_mpz_get_str(cuda_mpz_t *cuda_mpz, char *str, int bufsize) {
  int print_zeroes = 0; // don't print leading 0s
  int i;
  int str_index = 0;
  if (cuda_mpz_is_negative(cuda_mpz)) {
    str[0] = '-';
    str_index = 1;
  }

  #pragma unroll
  for (i = DIGITS_CAPACITY - 1; i >= 0; i--) {
    unsigned digit = cuda_mpz->digits[i];

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

__device__ __host__ inline void cuda_mpz_set(cuda_mpz_t *to, cuda_mpz_t *from) {
  unsigned i;

  #pragma unroll
  for (i = 0; i < from->words; i++) {// changes
    to->digits[i] = from->digits[i];
  }

  #pragma unroll
  for (; i < to->words; i++) {// changes
    to->digits[i] = 0;
  }

  to->sign = from->sign;
  to->words = from->words;
  to->bits = from->bits;
  //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;
}

__device__ __host__ inline void cuda_mpz_set_gmp(cuda_mpz_t *to, mpz_t from) {//changes
  int i;
  int word_count = 0;

  int src_size = from->_mp_size * 2;

  #pragma unroll
  for (i = src_size - 1; i >= 0; i--) {
	if((i & 1) == 0 ){
		to->digits[i] = (digit_t) (from->_mp_d[i/2] & 0xffffffff );
	}else{
		to->digits[i] = (digit_t) (from->_mp_d[i/2] >> 32);
	}
	if(to->digits[i] != 0 && word_count == 0){
		word_count = i + 1;
	}
  }

  #pragma unroll
  for (i = src_size; i < DIGITS_CAPACITY; i++) {
    to->digits[i] = 0;
  }

  to->sign = MPZ_NONNEGATIVE;//RSA only use positive numbers

  to->words = word_count;
  //finding the msb
  digit_t v = to->digits[word_count - 1];
  int msb = 0;

  while (v >>= 1) {
	  msb++;
  }

  to->bits = (word_count - 1) * LOG2_DIGIT_BASE + msb + 1;
  //to->words = (to->bits + LOG2_DIGIT_BASE - 1 ) / LOG2_DIGIT_BASE;
}

__host__ inline void cuda_mpz_set_str_host(cuda_mpz_t *cuda_mpz, const char *user_str) {//changes
  unsigned num_digits;
  unsigned i;
  int is_zero;

  #pragma unroll
  for (i = 0; i < DIGITS_CAPACITY; i++) cuda_mpz->digits[i] = 0;//changes

  const int bufsize = 1024;
  char buf[bufsize];
  memcpy(buf, user_str, strlen(user_str) + 1);
  buf[bufsize - 1] = (char) 0;
  char *str = &buf[0];

  /* Check if the provided number is negative */
  if (str[0] == '-') {
    cuda_mpz->sign = MPZ_NEGATIVE;
    str ++; // the number starts at the next character
  }
  else {
    cuda_mpz->sign = MPZ_NONNEGATIVE;
  }
  int len = strlen(str);
  int char_per_digit = LOG2_DIGIT_BASE / 4;
  num_digits = (len + char_per_digit - 1) / char_per_digit;
  CHECK_MEM(cuda_mpz, num_digits);

  digits_set_zero(cuda_mpz->digits);

  is_zero = true;
  #pragma unroll
  for (i = 0; i < num_digits; i ++) {
    str[len - i * char_per_digit] = (char) 0;
    char *start = str + (int) max(len - (i + 1) * char_per_digit, 0);
    digit_t d = strtol(start, NULL, 16);

    /* keep track of whether or not every digit is zero */
    is_zero = is_zero && (d == 0);

    /* parse the string backwards (little endian order) */
    cuda_mpz->digits[i] = d;
  }

  /* Just in case the user gives us -0 as input */
  if (is_zero) cuda_mpz->sign = MPZ_NONNEGATIVE;
}


__device__ __host__ inline void cuda_mpz_bitwise_truncate(cuda_mpz_t *dst, cuda_mpz_t *src) {//changes
  digit_t *src_digits = src->digits;
  digit_t *dst_digits = dst->digits;
  //unsigned capacity = src->capacity;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  //int ls_digits = capacity - rs_digits - 1;
  int ls_remainder = LOG2_DIGIT_BASE - ( RL & MOD_LOG2_DIGIT_BASE );
  //int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  dst->sign = src->sign;

  #pragma unroll
  for(int d_index = DIGITS_CAPACITY - 1; d_index > rs_digits; d_index--) {//constant time for specific rl
	  dst_digits[d_index] = 0;
  }

  dst_digits[rs_digits] = src_digits[rs_digits] & ( 0xffffffff >> ls_remainder);

  //pragma unroll
  for(int d_index = rs_digits - 1; d_index >= 0; d_index--) {//constant time for specific rl
	  dst_digits[d_index] = src_digits[d_index];
  }
}

__device__ __host__ inline void cuda_mpz_bitwise_truncate_eq(cuda_mpz_t *cuda_mpz) {//changes
  digit_t *digits = cuda_mpz->digits;
  //unsigned capacity = cuda_mpz->capacity;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  //int ls_digits = capacity - rs_digits - 1;
  int ls_remainder = LOG2_DIGIT_BASE - ( RL & MOD_LOG2_DIGIT_BASE );
  //int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;


  #pragma unroll
  for(int d_index = DIGITS_CAPACITY - 1; d_index > rs_digits; d_index--) {//constant time for specific rl
	digits[d_index] = 0;
  }

  digits[rs_digits] = digits[rs_digits] & ( 0xffffffff >> ls_remainder);


}


/**
 * @brief Perform dst := op1 * op2.
 *
 * @warning Assumes that all cuda_mpz_t parameters have been initialized.
 * @warning Assumes dst != op1 != op2
 */
__device__ __host__ inline void cuda_mpz_mult(cuda_mpz_t *dst, cuda_mpz_t *op1, cuda_mpz_t *op2) {
  unsigned op1_digit_count = digits_count(op1->digits);
  unsigned op2_digit_count = digits_count(op2->digits);
  unsigned capacity = max(op1_digit_count, op2_digit_count);

//  /* In multiplication, if the operand with the most digits has D digits,
//   * then the result of the addition will have at most 2D digits. */
//  CHECK_MEM(dst, 2*capacity);
//  CHECK_MEM(op1,   capacity);
//  CHECK_MEM(op2,   capacity); // changes

  /* Done by long_multiplication */
  /* digits_set_zero(dst->digits); */

  /* We pass in capacity as the number of digits rather that the actual
   * number of digits in each cuda_mpz_t struct. This is done because the
   * multiplication code has some assumptions and optimizations (e.g.
   * op1 and op2 to have the same number of digits) */
  digits_mult(dst->digits, op1->digits, op2->digits, capacity, DIGITS_CAPACITY);

  /* Compute the sign of the product */
  dst->sign = (op1->sign == op2->sign) ? MPZ_NONNEGATIVE : MPZ_NEGATIVE;
  if (MPZ_NEGATIVE == dst->sign && digits_is_zero(dst->digits, DIGITS_CAPACITY)) {
    dst->sign = MPZ_NONNEGATIVE;
  }

//  CHECK_SIGN(op1);
//  CHECK_SIGN(op2);
//  CHECK_SIGN(dst);

  #pragma unroll
  for (i = 0; i < num_digits; i++) {
    #pragma unroll
    for (j = 0; j < num_digits; j++) {
      unsigned k = i + j;
      digit_t carry = 0;
      digit_t prod;

      prod = mult(op2[i], op1[j], &carry);

      digits_add_across(product + k,     2*num_digits - k,     prod);
      digits_add_across(product + k + 1, 2*num_digits - k - 1, carry);
    }
  }
}

__device__ __host__ inline void clip(unsigned long long value,
                                     digit_t* result, digit_t *carry) {
  *carry  = (digit_t) (value / DIGIT_BASE); //FIXME
  *result = (digit_t) (value % DIGIT_BASE); //FIXME
  //printf("clip(%llu): result = %u, carry = %u\n", value, *result, *carry);
}

/**
 * @return
 *      < 0  if a < b
 *      = 0  if a = b
 *      > 0  if a > b
 *
 * @warning This function does not give any indication about the distance
 * between a and b, just the relative distance (<, >, =).
 */
#define MPZ_LESS    -1
#define MPZ_GREATER  1
#define MPZ_EQUAL    0
__device__ __host__ inline int cuda_mpz_compare(cuda_mpz_t *a, cuda_mpz_t *b) {
  int cmp;
  int negative;

  if (MPZ_NEGATIVE == a->sign && MPZ_NONNEGATIVE == b->sign) return MPZ_LESS;
  if (MPZ_NEGATIVE == b->sign && MPZ_NONNEGATIVE == a->sign) return MPZ_GREATER;

  /* At this point we know they have the same sign */
  cmp = digits_compare(a->digits, DIGITS_CAPACITY, b->digits, DIGITS_CAPACITY);
  negative = cuda_mpz_is_negative(a);

  if (cmp == 0) return MPZ_EQUAL;

  if (negative) {
    return (cmp > 0) ? MPZ_LESS : MPZ_GREATER;
  }
  else {
    return (cmp < 0) ? MPZ_LESS : MPZ_GREATER;
  }
}

/** @brief Return true if a >= b */
__device__ __host__ inline int cuda_mpz_gte(cuda_mpz_t *a, cuda_mpz_t *b) {
  return (cuda_mpz_compare(a, b) >= 0);
}

__device__ inline void cuda_mpz_print_str_device(cuda_mpz_t *cuda_mpz) {//changes
  int print_zeroes = 0; // don't print leading 0s

  if (cuda_mpz_is_negative(cuda_mpz)) {
	  printf("-");
  }

  #pragma unroll
  for (int i = DIGITS_CAPACITY - 1; i >= 0; i--) {
    unsigned digit = cuda_mpz->digits[i];

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

__device__ __host__ inline void cuda_mpz_bitwise_rshift_eq(cuda_mpz_t *cuda_mpz) {//changes
  digit_t *digits = cuda_mpz->digits;
  //unsigned capacity = cuda_mpz->capacity;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  int rs_remainder = RL & MOD_LOG2_DIGIT_BASE;

  #pragma unroll
  for(int d_index = 0; d_index < DIGITS_CAPACITY - 1 - rs_digits; d_index++) {//constant time for specific rl
	digits[d_index] = ( digits[d_index + rs_digits] >> rs_remainder ) | ( digits[d_index + rs_digits + 1] << ( LOG2_DIGIT_BASE - rs_remainder ) );
  }

  digits[DIGITS_CAPACITY - 1 - rs_digits] = digits[DIGITS_CAPACITY - 1] >> rs_remainder;

  #pragma unroll
  for(int d_index = DIGITS_CAPACITY - rs_digits; d_index <= DIGITS_CAPACITY - 1; d_index++) {//constant time for specific rl
  	digits[d_index] = 0;
  }
}

__device__ __host__ inline void cuda_mpz_bitwise_rshift(cuda_mpz_t *dst, cuda_mpz_t *src) {//changes
  digit_t *src_digits = src->digits;
  digit_t *dst_digits = dst->digits;
  //unsigned capacity = src->capacity;

  int rs_digits = RL >> LOG2_LOG2_DIGIT_BASE;
  int rs_remainder = RL & MOD_LOG2_DIGIT_BASE;

  dst->sign = src->sign;

  #pragma unroll
  for(int d_index = 0; d_index < DIGITS_CAPACITY - 1 - rs_digits; d_index++) {//constant time for specific rl
	  dst_digits[d_index] = ( src_digits[d_index + rs_digits] >> rs_remainder ) | ( src_digits[d_index + rs_digits + 1] << ( LOG2_DIGIT_BASE - rs_remainder ) );
  }

  dst_digits[DIGITS_CAPACITY - 1 - rs_digits] = src_digits[DIGITS_CAPACITY - 1] >> rs_remainder;

  #pragma unroll
  for(int d_index = DIGITS_CAPACITY - rs_digits; d_index <= DIGITS_CAPACITY - 1; d_index++) {//constant time for specific rl
	  dst_digits[d_index] = 0;
  }
}

__device__ __host__ inline digit_t cuda_mpz_get_last_digit(cuda_mpz_t *cuda_mpz) {//changes
	return cuda_mpz->digits[0];
}


#endif /* __418_MPZ_H__ */
