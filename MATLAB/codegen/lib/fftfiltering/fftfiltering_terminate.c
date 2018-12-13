/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fftfiltering_terminate.c
 *
 * Code generation for function 'fftfiltering_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fftfiltering.h"
#include "fftfiltering_terminate.h"
#include "fftfiltering_data.h"

/* Function Definitions */
void fftfiltering_terminate(void)
{
  omp_destroy_nest_lock(&emlrtNestLockGlobal);
}

/* End of code generation (fftfiltering_terminate.c) */
