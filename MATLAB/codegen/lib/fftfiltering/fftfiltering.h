/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fftfiltering.h
 *
 * Code generation for function 'fftfiltering'
 *
 */

#ifndef FFTFILTERING_H
#define FFTFILTERING_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "fftfiltering_types.h"

/* Function Declarations */
extern void fftfiltering(const emxArray_real32_T *b, emxArray_real32_T *x, float
  nfft, emxArray_creal32_T *y);

#endif

/* End of code generation (fftfiltering.h) */
