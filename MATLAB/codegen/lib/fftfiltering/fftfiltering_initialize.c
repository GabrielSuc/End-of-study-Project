/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fftfiltering_initialize.c
 *
 * Code generation for function 'fftfiltering_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fftfiltering.h"
#include "fftfiltering_initialize.h"
#include "fftfiltering_data.h"

/* Function Definitions */
void fftfiltering_initialize(void)
{
  rt_InitInfAndNaN(8U);
  omp_init_nest_lock(&emlrtNestLockGlobal);
}

/* End of code generation (fftfiltering_initialize.c) */
