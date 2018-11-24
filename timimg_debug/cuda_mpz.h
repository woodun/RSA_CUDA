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
  char     sign;
} cuda_mpz_t;

__host__ inline char* cuda_mpz_get_str(cuda_mpz_t *cuda_mpz, char *str, int bufsize);

/**
 * @brief Check that the cuda_mpz_t struct has enough memory to store __capacity
 * digits.
 */
#ifndef __CUDACC__
#define CHECK_MEM(__cuda_mpz, __capacity) \
  do {                                                                  \
    if ((__cuda_mpz)->capacity < (__capacity)) {                             \
      printf("MPZ memory error at %s:%d.\n", __func__, __LINE__);       \
      printf("\tcuda_mpz capacity: %u, requested capacity %u\n",             \
             (__cuda_mpz)->capacity, (__capacity));                          \
      assert(0);                                                        \
    }                                                                   \
  } while (0)
#else
#define CHECK_MEM(__cuda_mpz, __capacity)
#endif

/**
 * @brief Sanity check the sign.
 */
#ifndef __CUDACC__
#define CHECK_SIGN(__cuda_mpz) \
  do {                                                                  \
    if (digits_is_zero((__cuda_mpz)->digits, (__cuda_mpz)->capacity) &&           \
        (__cuda_mpz)->sign != MPZ_NONNEGATIVE) {                             \
      printf("MPZ Sign Error at %s:%d: Value is 0 but sign is %d.\n",   \
             __func__, __LINE__, (__cuda_mpz)->sign);                        \
      assert(0);\
    }                                                                   \
    if ((__cuda_mpz)->sign != MPZ_NEGATIVE &&                                \
        (__cuda_mpz)->sign != MPZ_NONNEGATIVE) {                             \
      printf("MPZ Sign Error at %s:%d: Invalid sign %d.\n",             \
             __func__, __LINE__, (__cuda_mpz)->sign);                        \
      assert(0);\
    }                                                                   \
  } while (0)
#else
#define CHECK_SIGN(__cuda_mpz)
#endif

/**
 * @brief Do some sanity checking on the cuda_mpz_t sign field.
 */
#ifndef __CUDACC__
#define CHECK_STRS(s1, s2)                                            \
  do {                                                                \
    if (strcmp(s1, s2)) {                                             \
      printf("Input string %s became %s!\n", s1, s2);                 \
    }                                                                 \
  } while (0)
#else
#define CHECK_STRS(s1, s2)
#endif


__device__ __host__ inline int cuda_mpz_is_negative(cuda_mpz_t *cuda_mpz) {
  return (cuda_mpz->sign == MPZ_NEGATIVE);
}

__device__ __host__ inline void cuda_mpz_negate(cuda_mpz_t *cuda_mpz) {
  cuda_mpz->sign = (cuda_mpz->sign == MPZ_NEGATIVE) ? MPZ_NONNEGATIVE : MPZ_NEGATIVE;
  if (digits_is_zero(cuda_mpz->digits, cuda_mpz->capacity)) {
    cuda_mpz->sign = MPZ_NONNEGATIVE;
  }
}

/**
 * @brief Initialize an cuda_mpz struct to 0.
 */
__device__ __host__ inline void cuda_mpz_init(cuda_mpz_t *cuda_mpz) {
  cuda_mpz->capacity = DIGITS_CAPACITY;
  /* digits_set_zero(cuda_mpz->digits); */
  cuda_mpz->sign = MPZ_NONNEGATIVE; //changes
}

/**
 * @brief Assign an cuda_mpz_t struct to the value of another cuda_mpz_t struct.
 */
__device__ __host__ inline void cuda_mpz_set(cuda_mpz_t *to, cuda_mpz_t *from) {
  unsigned i;

  #pragma unroll
  for (i = 0; i < DIGITS_CAPACITY; i++) {// changes
    digit_t d = (i < DIGITS_CAPACITY) ? from->digits[i] : 0;//changes
    to->digits[i] = d;
  }

  to->sign = from->sign;

  //CHECK_SIGN(to);
  //CHECK_SIGN(from);
}

__device__ __host__ inline void cuda_mpz_set_gmp(cuda_mpz_t *to, mpz_t from) {//changes
  unsigned i;

  /////////////////////////////////////////////////////todo: also copy length from from
  int src_size = from->_mp_size * 2;

  #pragma unroll
  for (i = 0; i < src_size; i++) {
	if((i & 1) == 0 ){
		to->digits[i] = (digit_t) (from->_mp_d[i/2] & 0xffffffff );
	}else{
		to->digits[i] = (digit_t) (from->_mp_d[i/2] >> 32);
	}
  }

  #pragma unroll
  for (i = src_size; i < DIGITS_CAPACITY; i++) {
    to->digits[i] = 0;
  }

  to->sign = MPZ_NONNEGATIVE;//RSA only use positive numbers
}

/**
 * @brief Set the cuda_mpz integer to the provided integer.
 */
__device__ __host__ inline void cuda_mpz_set_i(cuda_mpz_t *cuda_mpz, int z) {
  cuda_mpz->sign = (z < 0) ? MPZ_NEGATIVE : MPZ_NONNEGATIVE;
  digits_set_lui(cuda_mpz->digits, abs(z));
}

/**
 * @brief Set the cuda_mpz integer to the provided integer.
 */
__device__ __host__ inline void cuda_mpz_set_lui(cuda_mpz_t *cuda_mpz, unsigned long z) {
  cuda_mpz->sign = MPZ_NONNEGATIVE;
  digits_set_lui(cuda_mpz->digits, z);
}

__device__ __host__ inline void cuda_mpz_set_llui(cuda_mpz_t *cuda_mpz, unsigned long long z) {// changes
  cuda_mpz->sign = MPZ_NONNEGATIVE;
  digits_set_llui(cuda_mpz->digits, z);
}

/**
 * @brief Set the cuda_mpz integer to the provided integer.
 */
__device__ __host__ inline void cuda_mpz_set_ui(cuda_mpz_t *cuda_mpz, unsigned z) {
  cuda_mpz->sign = MPZ_NONNEGATIVE;
  digits_set_lui(cuda_mpz->digits, z);
}

/**
 * @brief Set the cuda_mpz integer based on the provided (hex) string.
 */
__device__ __host__ inline void cuda_mpz_set_str(cuda_mpz_t *cuda_mpz, const char *user_str) {
  unsigned num_digits;
  unsigned i;
  int is_zero;

  #pragma unroll
  for (i = 0; i < DIGITS_CAPACITY; i++) cuda_mpz->digits[i] = 0;

  const int bufsize = 1024;
  char buf[bufsize];
  memcpy(buf, user_str, bufsize);
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
  int len = cuda_strlen(str);
  int char_per_digit = LOG2_DIGIT_BASE / 4;
  num_digits = (len + char_per_digit - 1) / char_per_digit;
  CHECK_MEM(cuda_mpz, num_digits);

  digits_set_zero(cuda_mpz->digits);

  is_zero = true;
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

#if 0
  cuda_mpz_get_str(cuda_mpz, buf, bufsize);
  CHECK_STRS(user_str, buf);
#endif
}

/**
 * @brief Set the cuda_mpz integer based on the provided (hex) string.
 */
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

#if 0
  cuda_mpz_get_str(cuda_mpz, buf, bufsize);
  CHECK_STRS(user_str, buf);
#endif
}

__device__ __host__ inline void cuda_mpz_get_binary_str(cuda_mpz_t *cuda_mpz, char *str, unsigned s) {
  (void) cuda_mpz;
  (void) str;
  (void) s;
}

/**
 * @brief Destroy the cuda_mpz_t struct.
 *
 * @deprecated
 */
__device__ __host__ inline void cuda_mpz_destroy(cuda_mpz_t *cuda_mpz) {
  (void) cuda_mpz;
}

/**
 * @brief Add two multiple precision integers.
 *
 *      dst := op1 + op2
 *
 * @warning It is assumed that all cuda_mpz_t parameters have been initialized.
 * @warning Assumes dst != op1 != op2
 */
__device__ __host__ inline void cuda_mpz_add(cuda_mpz_t *dst, cuda_mpz_t *op1, cuda_mpz_t *op2) {
//#ifdef __CUDACC__
//  unsigned op1_digit_count = digits_count(op1->digits);
//  unsigned op2_digit_count = digits_count(op2->digits);
//
//  /* In addition, if the operand with the most digits has D digits, then
//   * the result of the addition will have at most D + 1 digits. */
//  unsigned capacity = max(op1_digit_count, op2_digit_count) + 1;
//
//  /* Make sure all of the cuda_mpz structs have enough memory to hold all of
//   * the digits. We will be doing 10's complement so everyone needs to
//   * have enough digits. */
//  CHECK_MEM(dst, capacity);
//  CHECK_MEM(op1, capacity);
//  CHECK_MEM(op2, capacity);
//#endif //changes

  digits_set_zero(dst->digits);

  /* If both are negative, treate them as positive and negate the result */
  if (cuda_mpz_is_negative(op1) && cuda_mpz_is_negative(op2)) {
    digits_add(dst->digits, DIGITS_CAPACITY,
               op1->digits, DIGITS_CAPACITY,
               op2->digits, DIGITS_CAPACITY);//changes ////todo: using real length here?
    dst->sign = MPZ_NEGATIVE;
  }
  /* one or neither are negative */
  else {
    digit_t carry_out;

    /* Perform 10's complement on negative numbers before adding */
    if (cuda_mpz_is_negative(op1)) digits_complement(op1->digits, DIGITS_CAPACITY);
    if (cuda_mpz_is_negative(op2)) digits_complement(op2->digits, DIGITS_CAPACITY); // changes

    carry_out = digits_add(dst->digits, DIGITS_CAPACITY,
                           op1->digits, DIGITS_CAPACITY,
                           op2->digits, DIGITS_CAPACITY); // changes

    /* If there is no carryout, the result is negative */
    if (carry_out == 0 && (cuda_mpz_is_negative(op1) || cuda_mpz_is_negative(op2))) {
      digits_complement(dst->digits, DIGITS_CAPACITY);// changes
      dst->sign = MPZ_NEGATIVE;
    }
    /* Otherwise, the result is non-negative */
    else {
      dst->sign = MPZ_NONNEGATIVE;
    }

    /* Undo the 10s complement after adding */
    if (cuda_mpz_is_negative(op1)) digits_complement(op1->digits, DIGITS_CAPACITY);
    if (cuda_mpz_is_negative(op2)) digits_complement(op2->digits, DIGITS_CAPACITY);
  }

//  CHECK_SIGN(op1);
//  CHECK_SIGN(op2);
//  CHECK_SIGN(dst);
}

__device__ __host__ inline void cuda_mpz_bitwise_and(cuda_mpz_t *dst, cuda_mpz_t *op1, cuda_mpz_t *op2) {//changes
//#ifdef __CUDACC__
//  unsigned op1_digit_count = digits_count(op1->digits);
//  unsigned op2_digit_count = digits_count(op2->digits);
//
//  unsigned capacity = min(op1_digit_count, op2_digit_count);
//#endif

  digits_set_zero(dst->digits);

  /*
  if (cuda_mpz_is_negative(op1) || cuda_mpz_is_negative(op2)) {//do not accept negative numbers, will not happen in RSA.
	 //assert(0);
	 //return;
  }
  */

  digits_bitwise_and(dst->digits, DIGITS_CAPACITY, op1->digits, DIGITS_CAPACITY, op2->digits, DIGITS_CAPACITY);
}

__device__ __host__ inline void cuda_mpz_bitwise_truncate(cuda_mpz_t *dst, cuda_mpz_t *src, int n_bits) {//changes
  digit_t *src_digits = src->digits;
  digit_t *dst_digits = dst->digits;
  //unsigned capacity = src->capacity;

  int rs_digits = n_bits >> LOG2_LOG2_DIGIT_BASE;
  //int ls_digits = capacity - rs_digits - 1;
  int ls_remainder = LOG2_DIGIT_BASE - ( n_bits & MOD_LOG2_DIGIT_BASE );
  //int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  /*
  if(rs_digits >= capacity){//this will not happen in the rsa, but in the general case should be added.
	  return;
  }
  if(n_bits <= 0){//this will not happen in the rsa, but in the general case should be added.
	  digits_set_zero(digits);
	  return;
  }
  */

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

__device__ __host__ inline void cuda_mpz_bitwise_truncate_eq(cuda_mpz_t *cuda_mpz, int n_bits) {//changes
  digit_t *digits = cuda_mpz->digits;
  //unsigned capacity = cuda_mpz->capacity;

  int rs_digits = n_bits >> LOG2_LOG2_DIGIT_BASE;
  //int ls_digits = capacity - rs_digits - 1;
  int ls_remainder = LOG2_DIGIT_BASE - ( n_bits & MOD_LOG2_DIGIT_BASE );
  //int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  /*
  if(rs_digits >= capacity){//this will not happen in the rsa, but in the general case should be added.
	  return;
  }
  if(n_bits <= 0){//this will not happen in the rsa, but in the general case should be added.
	  digits_set_zero(digits);
	  return;
  }
  */

  #pragma unroll
  for(int d_index = DIGITS_CAPACITY - 1; d_index > rs_digits; d_index--) {//constant time for specific rl
	digits[d_index] = 0;
  }

  digits[rs_digits] = digits[rs_digits] & ( 0xffffffff >> ls_remainder);
}

__device__ __host__ inline void cuda_mpz_addeq(cuda_mpz_t *op1, cuda_mpz_t *op2) {

  /* If both are negative, treate them as positive and negate the result */
  if (cuda_mpz_is_negative(op1) && cuda_mpz_is_negative(op2)) {
    digits_addeq(op1->digits,DIGITS_CAPACITY,
               op2->digits, DIGITS_CAPACITY);
    op1->sign = MPZ_NEGATIVE;
  }
  /* one or neither are negative */
  else {
    digit_t carry_out;

    /* Perform 10's complement on negative numbers before adding */
    if (cuda_mpz_is_negative(op1)) digits_complement(op1->digits, DIGITS_CAPACITY);
    if (cuda_mpz_is_negative(op2)) digits_complement(op2->digits, DIGITS_CAPACITY);

    carry_out = digits_addeq(op1->digits, DIGITS_CAPACITY,
                             op2->digits, DIGITS_CAPACITY);

    /* If there is no carryout, the result is negative */
    if (carry_out == 0 && (cuda_mpz_is_negative(op1) || cuda_mpz_is_negative(op2))) {
      digits_complement(op1->digits, DIGITS_CAPACITY);
      op1->sign = MPZ_NEGATIVE;
    }
    /* Otherwise, the result is non-negative */
    else {
      op1->sign = MPZ_NONNEGATIVE;
    }

    /* Undo the 10s complement after adding */
    if (cuda_mpz_is_negative(op2)) digits_complement(op2->digits, DIGITS_CAPACITY);
  }

//  CHECK_SIGN(op1);
//  CHECK_SIGN(op2);
}

/**
 * @brief Perform dst := op1 - op2.
 *
 * @warning Assumes that all cuda_mpz_t parameters have been initialized.
 * @warning Assumes dst != op1 != op2
 */
__device__ __host__ inline void cuda_mpz_sub(cuda_mpz_t *dst, cuda_mpz_t *op1, cuda_mpz_t *op2) {
  cuda_mpz_negate(op2);
  cuda_mpz_add(dst, op1, op2);
  cuda_mpz_negate(op2);
}

/**
 * @brief Perform op1 -= op2.
 *
 * @warning Assumes that all cuda_mpz_t parameters have been initialized.
 * @warning Assumes op1 != op2
 */
__device__ __host__ inline void cuda_mpz_subeq(cuda_mpz_t *op1, cuda_mpz_t *op2) {
  cuda_mpz_negate(op2);
  cuda_mpz_addeq(op1, op2);
  cuda_mpz_negate(op2);
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

/** @brief Return true if a == b */
__device__ __host__ inline int cuda_mpz_equal(cuda_mpz_t *a, cuda_mpz_t *b) {
  return (cuda_mpz_compare(a, b) == 0);
}
/** @brief Return true if a == 1 */
__device__ __host__ inline int cuda_mpz_equal_one(cuda_mpz_t *a) {
  if (MPZ_NEGATIVE == a->sign) {
    return false;
  }
  return digits_equal_one(a->digits, DIGITS_CAPACITY);
}
/** @brief Return true if a < b */
__device__ __host__ inline int cuda_mpz_lt(cuda_mpz_t *a, cuda_mpz_t *b) {
  return (cuda_mpz_compare(a, b) < 0);
}
/** @brief Return true if a <= b */
__device__ __host__ inline int cuda_mpz_lte(cuda_mpz_t *a, cuda_mpz_t *b) {
  return (cuda_mpz_compare(a, b) <= 0);
}
/** @brief Return true if a > b */
__device__ __host__ inline int cuda_mpz_gt(cuda_mpz_t *a, cuda_mpz_t *b) {
  return (cuda_mpz_compare(a, b) > 0);
}
/** @brief Return true if a == 1 */
__device__ __host__ inline int cuda_mpz_gt_one(cuda_mpz_t *a) {
  if (MPZ_NEGATIVE == a->sign) {
    return false;
  }
  return digits_gt_one(a->digits, DIGITS_CAPACITY);
}
/** @brief Return true if a >= b */
__device__ __host__ inline int cuda_mpz_gte(cuda_mpz_t *a, cuda_mpz_t *b) {
  return (cuda_mpz_compare(a, b) >= 0);
}

/**
 * @brief Return the string representation of the integer represented by the
 * cuda_mpz_t struct.
 *
 * @warning If buf is NULL, the string is dynamically allocated and must
 * therefore be freed by the user.
 */
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

__host__ inline void cuda_mpz_print(cuda_mpz_t *cuda_mpz) {
  char str[1024];
  cuda_mpz_get_str(cuda_mpz, str, 1024);
  printf("%s\n", str);
}

__host__ __device__ inline void cuda_mpz_set_bit(cuda_mpz_t *cuda_mpz, unsigned bit_offset,
                                            unsigned bit) {//changes
  digits_set_bit(cuda_mpz->digits, bit_offset, bit);

  if (MPZ_NEGATIVE == cuda_mpz->sign && bit == 0 &&
      digits_is_zero(cuda_mpz->digits, DIGITS_CAPACITY)) {
    cuda_mpz->sign = MPZ_NONNEGATIVE;
  }
}

__device__ __host__ inline void cuda_mpz_bit_lshift(cuda_mpz_t *cuda_mpz) {
  bits_lshift(cuda_mpz->digits, DIGITS_CAPACITY);

  if (MPZ_NEGATIVE == cuda_mpz->sign && digits_is_zero(cuda_mpz->digits, DIGITS_CAPACITY)) {
    cuda_mpz->sign = MPZ_NONNEGATIVE;
  }
}

__device__ __host__ inline void cuda_mpz_bit_rshift(cuda_mpz_t *cuda_mpz, int n_bits) {
  #pragma unroll
  for(int i=0;i<n_bits;i++) bits_rshift(cuda_mpz->digits, DIGITS_CAPACITY);
}

__device__ __host__ inline void cuda_mpz_bitwise_rshift_eq(cuda_mpz_t *cuda_mpz, int n_bits) {//changes
  digit_t *digits = cuda_mpz->digits;
  //unsigned capacity = cuda_mpz->capacity;

  int rs_digits = n_bits >> LOG2_LOG2_DIGIT_BASE;
  int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  /*
  if(rs_digits >= capacity - 1){//this will not happen in the rsa, but in the general case should be added.
  	  digit_t first_digit = digits[capacity - 1];
	  digits_set_zero(digits);
	  digits[0] = first_digit >> (rs_digits - capacity + 1) * LOG2_DIGIT_BASE + rs_remainder;
	  return;
  }
  */

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

__device__ __host__ inline void cuda_mpz_bitwise_rshift(cuda_mpz_t *dst, cuda_mpz_t *src, int n_bits) {//changes
  digit_t *src_digits = src->digits;
  digit_t *dst_digits = dst->digits;
  //unsigned capacity = src->capacity;

  int rs_digits = n_bits >> LOG2_LOG2_DIGIT_BASE;
  int rs_remainder = n_bits & MOD_LOG2_DIGIT_BASE;

  /*
  if(rs_digits >= capacity - 1){//this will not happen in the rsa, but in the general case should be added.
  	  digit_t first_digit = digits[capacity - 1];
	  digits_set_zero(digits);
	  digits[0] = first_digit >> (rs_digits - capacity + 1) * LOG2_DIGIT_BASE + rs_remainder;
	  return;
  }
  */

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

__device__ __host__ inline int cuda_mpz_is_zero(cuda_mpz_t *cuda_mpz) {
  unsigned d_index = 0;

  #pragma unroll
  for (d_index = 0; d_index < DIGITS_CAPACITY; d_index++) {
    //printf("d: %d\n", cuda_mpz->digits[d_index]);
    if(cuda_mpz->digits[d_index] != 0) return 0;
  }
  return 1; //true
}

__device__ __host__ inline void cuda_mpz_div(cuda_mpz_t *q, cuda_mpz_t *r, cuda_mpz_t *n,
                                            cuda_mpz_t *d) {
  unsigned n_digit_count = digits_count(n->digits);
  unsigned num_bits;
  int i;
  int nsign = n->sign;
  int dsign = d->sign;

  num_bits = n_digit_count * LOG2_DIGIT_BASE;

  cuda_mpz_set_ui(q, 0);
  cuda_mpz_set_ui(r, 0);

  n->sign = MPZ_NONNEGATIVE;
  d->sign = MPZ_NONNEGATIVE;

  if (cuda_mpz_gt(n, d)) {

	#pragma unroll
    for (i = num_bits - 1; i >= 0; i--) {
      unsigned n_i;

      // r = r << 1
      cuda_mpz_bit_lshift(r);

      // r(0) = n(i)
      n_i = digits_bit_at(n->digits, i);
      cuda_mpz_set_bit(r, 0, n_i);

      // if (r >= d)
      if (cuda_mpz_gte(r, d)) {
        // r = r - d
        cuda_mpz_subeq(r, d);

        // q(i) = 1
        //printf("Setting bit %d of q to 1\n", i);
        //printf("\tBefore: "); cuda_mpz_print(q); printf("\n");
        cuda_mpz_set_bit(q, i, 1);
        //printf("\tAfter: "); cuda_mpz_print(q); printf("\n");
      }
    }

    /* Compute the sign of the division */
    q->sign = (nsign == dsign) ? MPZ_NONNEGATIVE : MPZ_NEGATIVE;
    if (MPZ_NEGATIVE == q->sign && digits_is_zero(q->digits, DIGITS_CAPACITY)) {
      q->sign = MPZ_NONNEGATIVE;
    }
  }
  else {
    // quotient = 0
    cuda_mpz_set_ui(q, 0);
    // remainder = numerator
    cuda_mpz_set(r, n);
  }

  n->sign = nsign;
  d->sign = dsign;

//  CHECK_SIGN(q);
//  CHECK_SIGN(r);
//  CHECK_SIGN(n);
//  CHECK_SIGN(d);
}

/**
 * @brief Compute the GCD of op1 and op2.
 *
 * Euclidean Algorithm:
 *
 *    while (b != 0) {
 *      t := b
 *      b := a % b
 *      a := t
 *    }
 *    gcd = a
 */
__device__ __inline__ void cuda_mpz_gcd_tmp(cuda_mpz_t *gcd, cuda_mpz_t *op1, cuda_mpz_t *op2,
                                       // tmps
                                       cuda_mpz_t *tmp1, cuda_mpz_t *tmp2,
                                       cuda_mpz_t *tmp3) {
  cuda_mpz_t *a = gcd;
  cuda_mpz_t *b = tmp1;
  cuda_mpz_t *mod = tmp2;
  cuda_mpz_t *quo = tmp3;

  int compare = cuda_mpz_compare(op1, op2);

  cuda_mpz_set(a, (compare > 0) ? op1 : op2);
  cuda_mpz_set(b, (compare > 0) ? op2 : op1);

  #pragma unroll
  while (!digits_is_zero(b->digits, DIGITS_CAPACITY)) {
    cuda_mpz_div(quo, mod, a, b);
    cuda_mpz_set(a, b);
    cuda_mpz_set(b, mod);
  }
}

__device__ __inline__ void cuda_mpz_gcd(cuda_mpz_t *gcd, cuda_mpz_t *op1, cuda_mpz_t *op2) {
  cuda_mpz_t tmp1;
  cuda_mpz_t tmp2;
  cuda_mpz_t tmp3;

  cuda_mpz_init(&tmp1);
  cuda_mpz_init(&tmp2);
  cuda_mpz_init(&tmp3);

  cuda_mpz_gcd_tmp(gcd, op1, op2, &tmp1, &tmp2, &tmp3);
}

/**
 * Using exponentiation by squaring algorithm:
 *
 *  function modular_pow(base, exponent, modulus)
 *    result := 1
 *    while exponent > 0
 *      if (exponent mod 2 == 1):
 *         result := (result * base) mod modulus
 *      exponent := exponent >> 1
 *      base = (base * base) mod modulus
 *    return result
 */
__device__ __inline__ void cuda_mpz_powmod_tmp(cuda_mpz_t *result, cuda_mpz_t *base,
                                          cuda_mpz_t *exp, cuda_mpz_t *mod,
                                          // temps
                                          cuda_mpz_t *tmp1, cuda_mpz_t *tmp2,
                                          cuda_mpz_t *tmp3) {
  unsigned iteration;

  cuda_mpz_t *b = tmp3;

  // result = 1
  cuda_mpz_set_ui(result, 1);

  // _base = base % mod
  cuda_mpz_set(tmp1, base);
  cuda_mpz_div(tmp2, b, tmp1, mod);

  iteration = 0;
  #pragma unroll
  while (!bits_is_zero(exp->digits, DIGITS_CAPACITY, iteration)) {
    // if (binary_exp is odd)
    if (digits_bit_at(exp->digits, iteration) == 1) {
      // result = (result * base) % mod
      cuda_mpz_mult(tmp1, result, b);
      cuda_mpz_div(tmp2, result, tmp1, mod);
    }

    // binary_exp = binary_exp >> 1
    iteration++;

    // base = (base * base) % mod
    cuda_mpz_set(tmp1, b);
    cuda_mpz_mult(tmp2, b, tmp1);
    cuda_mpz_div(tmp1, b, tmp2, mod);
  }
}

__device__ __inline__ void cuda_mpz_powmod(cuda_mpz_t *result, cuda_mpz_t *base,
                                      cuda_mpz_t *exp, cuda_mpz_t *mod) {
  cuda_mpz_t tmp1;
  cuda_mpz_t tmp2;
  cuda_mpz_t tmp3;

  cuda_mpz_init(&tmp1);
  cuda_mpz_init(&tmp2);
  cuda_mpz_init(&tmp3);

  cuda_mpz_powmod_tmp(result, base, exp, mod, &tmp1, &tmp2, &tmp3);
}



__device__ __inline__ void cuda_mpz_pow(cuda_mpz_t *result, cuda_mpz_t *base, unsigned exponent) {
  cuda_mpz_t tmp;
  unsigned int i;

  cuda_mpz_init(&tmp);

  // result = 1
  cuda_mpz_set_ui(result, 1);
  #pragma unroll
  for (i = 0; i < exponent; i++) {
    // result *= base
    cuda_mpz_mult(&tmp, result, base);
    cuda_mpz_set(result, &tmp);
  }
}

/**
 * @brief Compute a += i
 */
__device__ __inline__ void cuda_mpz_addeq_i(cuda_mpz_t *a, int i) {
  if (0 == i) return;

  if (i < 0) {
    digits_complement(a->digits, DIGITS_CAPACITY);
    digits_add_across(a->digits, DIGITS_CAPACITY, -i);
    digits_complement(a->digits, DIGITS_CAPACITY);
  }
  else {
    digits_add_across(a->digits, DIGITS_CAPACITY, i);
  }
}

/**
 * @brief Compute result = a * i
 */
__device__ __inline__ void cuda_mpz_mult_u(cuda_mpz_t *result, cuda_mpz_t *a, unsigned i) {
  digits_mult_u(result->digits, a->digits, i);

  result->sign = a->sign;
  if (0 == i) result->sign = MPZ_NONNEGATIVE;

//  CHECK_SIGN(result);
}

#endif /* __418_MPZ_H__ */
