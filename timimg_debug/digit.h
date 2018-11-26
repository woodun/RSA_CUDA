/**
 * @file digit.c
 *
 * @brief Library of functions that operate on arrays of digits.
 *
 * Arrays of digits are assumed to be in little endian order.
 *
 * @author David Matlack (dmatlack)
 */
#ifndef __418_DIGIT_H__
#define __418_DIGIT_H__

#include "compile.h"

#define LOG2_DIGIT_BASE     32
#define DIGIT_BASE          ((unsigned long long) 1 << (LOG2_DIGIT_BASE))
#define DIGITS_CAPACITY     8 //changes: make enough space for large input
#define MOD_DIGIT_BASE      0xffffffff//changes
#define MOD_LOG2_DIGIT_BASE     31 //changes
#define LOG2_LOG2_DIGIT_BASE 5 //changes
#define RL 70

typedef unsigned digit_t;

__device__ __host__ inline void clip(unsigned long long value, digit_t* result, digit_t *carry) {
  *carry  = (digit_t) (value >> LOG2_DIGIT_BASE);
  *result = (digit_t) (value & MOD_DIGIT_BASE);
}

__device__ __host__ inline digit_t mult(digit_t a, digit_t b, digit_t *carry) {

  unsigned long long tmp = ((unsigned long long) a) * ((unsigned long long) b) +
                           ((unsigned long long) *carry);
  digit_t result;

  clip(tmp, &result, carry);

  return result;
}

__device__ __host__ inline digit_t add(digit_t a, digit_t b, digit_t *carry) {

  unsigned long long tmp = ((unsigned long long) a) +
                           ((unsigned long long) b) +
                           ((unsigned long long) *carry);
  digit_t result;

  clip(tmp, &result, carry);

  return result;
}

__device__ __host__ inline digit_t digits_add_across(digit_t *digits, unsigned num_digits, digit_t d) {
  digit_t carry = d;
  unsigned i = 0;

  while (carry != 0 && i < num_digits) {
    digits[i] = add(digits[i], 0, &carry);
    i++;
  }

  return carry;
}

/**
 * @brief Compute product = op1 * op2 using the Long Multiplication
 * (Grade School Multiplication) Aglorithm.
 *
 * @warning It is assumed that op1 and op2 contain num_digits each
 * and product has room for at least 2*num_digits.
 *
 * @return The carry out.
 */
__device__ __host__ inline void long_multiplication(digit_t *product,
                                                    digit_t *op1,
                                                    digit_t *op2,
                                                    unsigned num_digits,
                                                    unsigned prod_capacity) {
  int is_op1_zero = true;
  int is_op2_zero = true;
  unsigned i, j;

  /* zero out product */
  #pragma unroll
  for (i = 0; i < 2*num_digits; i++) {
    if (i < num_digits) {
      is_op1_zero = (is_op1_zero) && (op1[i] == 0);
      is_op2_zero = (is_op1_zero) && (op1[i] == 0);
    }
    product[i] = 0;
  }
  #pragma unroll
  for (; i < DIGITS_CAPACITY; i ++) {
    product[i] = 0;
  }

  /* if either of the operands are zero, then their product is zero */
  if (is_op1_zero || is_op2_zero) return;

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

/**
 * @brief Comare the two arrays of digits.
 *
 * @return  < 0 if d1 < d2
 *          = 0 if d1 == d2
 *          > 0 if d1 > d2
 *
 * @warning The return value does NOT give any indication about the relative
 * distance between the two numbers. It ONLY indicates <, >, or =.
 */
__device__ __host__ inline int digits_compare(digit_t *digits1, unsigned num_d1,
                                       digit_t *digits2, unsigned num_d2) {
  unsigned max_digits = (num_d1 > num_d2) ? num_d1 : num_d2;
  int i;

  /* Iterate backwards so that we look at the most significant digits first */
  for (i = max_digits - 1; i >= 0; i--) {
    digit_t d1 = ((unsigned) i < num_d1) ? digits1[i] : 0;
    digit_t d2 = ((unsigned) i < num_d2) ? digits2[i] : 0;

    if (d1 < d2) return -1;
    if (d1 > d2) return  1;
  }

  return 0;
}

#endif /* __418_DIGIT_H__ */
