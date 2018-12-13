/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fft1.h
 *
 * Code generation for function 'fft1'
 *
 */

#ifndef FFT1_H
#define FFT1_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "fftfiltering_types.h"

/* Function Declarations */
extern void b_fft(const emxArray_creal32_T *x, int n, emxArray_creal32_T *y);
extern void b_r2br_r2dit_trig(const emxArray_real32_T *x, int n1_unsigned, const
  emxArray_real32_T *costab, const emxArray_real32_T *sintab, emxArray_creal32_T
  *y);
extern void dobluesteinfft(const emxArray_real32_T *x, int N2, int n1, const
  emxArray_real32_T *costab, const emxArray_real32_T *sintab, const
  emxArray_real32_T *sintabinv, emxArray_creal32_T *y);
extern void fft(const emxArray_real32_T *x, int n, emxArray_creal32_T *y);
extern void generate_twiddle_tables(int nRows, boolean_T useRadix2,
  emxArray_real32_T *costab, emxArray_real32_T *sintab, emxArray_real32_T
  *sintabinv);
extern void get_algo_sizes(int n1, boolean_T useRadix2, int *N2blue, int *nRows);

#endif

/* End of code generation (fft1.h) */
