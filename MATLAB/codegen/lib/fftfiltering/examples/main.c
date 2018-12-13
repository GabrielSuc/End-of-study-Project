/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "fftfiltering.h"
#include "main.h"
#include "fftfiltering_terminate.h"
#include "fftfiltering_emxAPI.h"
#include "fftfiltering_initialize.h"

/* Function Declarations */
static emxArray_real32_T *argInit_1xUnbounded_real32_T(void);
static float argInit_real32_T(void);
static emxArray_real32_T *c_argInit_UnboundedxUnbounded_r(void);
static void main_fftfiltering(void);

/* Function Definitions */
static emxArray_real32_T *argInit_1xUnbounded_real32_T(void)
{
  emxArray_real32_T *result;
  static int iv1[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real32_T(2, iv1);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[result->size[0] * idx1] = argInit_real32_T();
  }

  return result;
}

static float argInit_real32_T(void)
{
  return 0.0F;
}

static emxArray_real32_T *c_argInit_UnboundedxUnbounded_r(void)
{
  emxArray_real32_T *result;
  static int iv2[2] = { 2, 2 };

  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real32_T(2, iv2);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_real32_T();
    }
  }

  return result;
}

static void main_fftfiltering(void)
{
  emxArray_creal32_T *y;
  emxArray_real32_T *b;
  emxArray_real32_T *x;
  emxInitArray_creal32_T(&y, 2);

  /* Initialize function 'fftfiltering' input arguments. */
  /* Initialize function input argument 'b'. */
  b = argInit_1xUnbounded_real32_T();

  /* Initialize function input argument 'x'. */
  x = c_argInit_UnboundedxUnbounded_r();

  /* Call the entry-point 'fftfiltering'. */
  fftfiltering(b, x, argInit_real32_T(), y);
  emxDestroyArray_creal32_T(y);
  emxDestroyArray_real32_T(x);
  emxDestroyArray_real32_T(b);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  fftfiltering_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_fftfiltering();

  /* Terminate the application.
     You do not need to do this more than one time. */
  fftfiltering_terminate();
  return 0;
}

/* End of code generation (main.c) */
