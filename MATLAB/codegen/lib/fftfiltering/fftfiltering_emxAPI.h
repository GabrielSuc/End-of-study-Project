/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fftfiltering_emxAPI.h
 *
 * Code generation for function 'fftfiltering_emxAPI'
 *
 */

#ifndef FFTFILTERING_EMXAPI_H
#define FFTFILTERING_EMXAPI_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "fftfiltering_types.h"

/* Function Declarations */
extern emxArray_creal32_T *emxCreateND_creal32_T(int numDimensions, int *size);
extern emxArray_real32_T *emxCreateND_real32_T(int numDimensions, int *size);
extern emxArray_creal32_T *emxCreateWrapperND_creal32_T(creal32_T *data, int
  numDimensions, int *size);
extern emxArray_real32_T *emxCreateWrapperND_real32_T(float *data, int
  numDimensions, int *size);
extern emxArray_creal32_T *emxCreateWrapper_creal32_T(creal32_T *data, int rows,
  int cols);
extern emxArray_real32_T *emxCreateWrapper_real32_T(float *data, int rows, int
  cols);
extern emxArray_creal32_T *emxCreate_creal32_T(int rows, int cols);
extern emxArray_real32_T *emxCreate_real32_T(int rows, int cols);
extern void emxDestroyArray_creal32_T(emxArray_creal32_T *emxArray);
extern void emxDestroyArray_real32_T(emxArray_real32_T *emxArray);
extern void emxInitArray_creal32_T(emxArray_creal32_T **pEmxArray, int
  numDimensions);
extern void emxInitArray_real32_T(emxArray_real32_T **pEmxArray, int
  numDimensions);

#endif

/* End of code generation (fftfiltering_emxAPI.h) */
