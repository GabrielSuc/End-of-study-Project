/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fftfiltering_emxutil.h
 *
 * Code generation for function 'fftfiltering_emxutil'
 *
 */

#ifndef FFTFILTERING_EMXUTIL_H
#define FFTFILTERING_EMXUTIL_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "fftfiltering_types.h"

/* Function Declarations */
extern void emxEnsureCapacity_creal32_T(emxArray_creal32_T *emxArray, int
  oldNumel);
extern void emxEnsureCapacity_creal32_T1(emxArray_creal32_T *emxArray, int
  oldNumel);
extern void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real32_T(emxArray_real32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real32_T1(emxArray_real32_T *emxArray, int
  oldNumel);
extern void emxFree_creal32_T(emxArray_creal32_T **pEmxArray);
extern void emxFree_int8_T(emxArray_int8_T **pEmxArray);
extern void emxFree_real32_T(emxArray_real32_T **pEmxArray);
extern void emxInit_creal32_T(emxArray_creal32_T **pEmxArray, int numDimensions);
extern void emxInit_creal32_T1(emxArray_creal32_T **pEmxArray, int numDimensions);
extern void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions);
extern void emxInit_real32_T(emxArray_real32_T **pEmxArray, int numDimensions);
extern void emxInit_real32_T1(emxArray_real32_T **pEmxArray, int numDimensions);

#endif

/* End of code generation (fftfiltering_emxutil.h) */
