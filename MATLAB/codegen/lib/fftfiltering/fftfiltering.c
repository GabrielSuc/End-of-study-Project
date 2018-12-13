/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fftfiltering.c
 *
 * Code generation for function 'fftfiltering'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "fftfiltering.h"
#include "fftfiltering_emxutil.h"
#include "fft1.h"
#include "any.h"

/* Function Declarations */
static float rt_powf_snf(float u0, float u1);

/* Function Definitions */
static float rt_powf_snf(float u0, float u1)
{
  float y;
  float f0;
  float f1;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else {
    f0 = (float)fabs(u0);
    f1 = (float)fabs(u1);
    if (rtIsInfF(u1)) {
      if (f0 == 1.0F) {
        y = 1.0F;
      } else if (f0 > 1.0F) {
        if (u1 > 0.0F) {
          y = ((real32_T)rtInf);
        } else {
          y = 0.0F;
        }
      } else if (u1 > 0.0F) {
        y = 0.0F;
      } else {
        y = ((real32_T)rtInf);
      }
    } else if (f1 == 0.0F) {
      y = 1.0F;
    } else if (f1 == 1.0F) {
      if (u1 > 0.0F) {
        y = u0;
      } else {
        y = 1.0F / u0;
      }
    } else if (u1 == 2.0F) {
      y = u0 * u0;
    } else if ((u1 == 0.5F) && (u0 >= 0.0F)) {
      y = (float)sqrt(u0);
    } else if ((u0 < 0.0F) && (u1 > (float)floor(u1))) {
      y = ((real32_T)rtNaN);
    } else {
      y = (float)pow(u0, u1);
    }
  }

  return y;
}

void fftfiltering(const emxArray_real32_T *b, emxArray_real32_T *x, float nfft,
                  emxArray_creal32_T *y)
{
  int fcnOutput;
  int b_x;
  int b_fcnOutput;
  int i0;
  emxArray_creal32_T *B;
  emxArray_creal32_T *r0;
  int b_b[1];
  emxArray_real32_T c_b;
  int loop_ub;
  int iv0[1];
  emxArray_real32_T *c_x;
  int unnamed_idx_1;
  float out[2];
  int i1;
  float istart;
  emxArray_creal32_T *X;
  emxArray_creal32_T *Y;
  emxArray_int8_T *r1;
  emxArray_real32_T *d_x;
  emxArray_real32_T *costab;
  emxArray_real32_T *sintab;
  emxArray_real32_T *sintabinv;
  emxArray_creal32_T *b_y;
  float y_re;
  float u1;
  emxArray_real32_T *d_b;
  int i2;
  boolean_T useRadix2;
  int i3;
  int b_Y[1];
  emxArray_creal32_T c_Y;
  float B_re;
  float B_im;
  fcnOutput = x->size[0];
  if (x->size[0] == 1.0F) {
    b_x = x->size[1];
    i0 = x->size[0] * x->size[1];
    x->size[0] = b_x;
    x->size[1] = 1;
    emxEnsureCapacity_real32_T(x, i0);

    /*  turn row into a column */
  }

  b_fcnOutput = x->size[0];

  /*  nfft is given */
  /*  Cast to enforce precision rules */
  if (nfft < 1.0F) {
    nfft = 1.0F;
  }

  emxInit_creal32_T(&B, 2);
  emxInit_creal32_T1(&r0, 1);
  nfft = rt_powf_snf(2.0F, (float)ceil((float)log(nfft) / 0.693147182F));

  /*  force this to a power of 2 for speed */
  b_b[0] = b->size[1];
  c_b = *b;
  c_b.size = (int *)&b_b;
  c_b.numDimensions = 1;
  fft(&c_b, (int)nfft, r0);
  i0 = B->size[0] * B->size[1];
  B->size[0] = 1;
  B->size[1] = (int)nfft;
  emxEnsureCapacity_creal32_T(B, i0);
  loop_ub = (int)nfft;
  for (i0 = 0; i0 < loop_ub; i0++) {
    B->data[i0] = r0->data[i0];
  }

  /*  if length(b)==1 */
  /*       B = B(:);  % make sure fft of B is a column (might be a row if b is scalar) */
  /*  end */
  if (b->size[1] == 1.0F) {
    iv0[0] = 1;
    c_b = *b;
    c_b.size = (int *)&iv0;
    c_b.numDimensions = 1;
    fft(&c_b, (int)nfft, r0);
    unnamed_idx_1 = (int)(float)x->size[1];
    i0 = B->size[0] * B->size[1];
    B->size[0] = 1;
    B->size[1] = unnamed_idx_1;
    emxEnsureCapacity_creal32_T(B, i0);
    for (i0 = 0; i0 < unnamed_idx_1; i0++) {
      B->data[B->size[0] * i0] = r0->data[0];
    }

    /*  replicate the column B  */
  }

  if (x->size[1] == 1.0F) {
    emxInit_real32_T(&c_x, 2);
    b_x = x->size[0];
    unnamed_idx_1 = (int)(float)b->size[1];
    i0 = c_x->size[0] * c_x->size[1];
    c_x->size[0] = b_x;
    c_x->size[1] = unnamed_idx_1;
    emxEnsureCapacity_real32_T(c_x, i0);
    for (i0 = 0; i0 < unnamed_idx_1; i0++) {
      for (i1 = 0; i1 < b_x; i1++) {
        c_x->data[i1 + c_x->size[0] * i0] = x->data[i1];
      }
    }

    i0 = x->size[0] * x->size[1];
    x->size[0] = c_x->size[0];
    x->size[1] = c_x->size[1];
    emxEnsureCapacity_real32_T(x, i0);
    loop_ub = c_x->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_x = c_x->size[0];
      for (i1 = 0; i1 < b_x; i1++) {
        x->data[i1 + x->size[0] * i0] = c_x->data[i1 + c_x->size[0] * i0];
      }
    }

    emxFree_real32_T(&c_x);

    /*  replicate the column x  */
  }

  for (i0 = 0; i0 < 2; i0++) {
    out[i0] = (float)x->size[i0];
  }

  i0 = y->size[0] * y->size[1];
  y->size[0] = (int)out[0];
  y->size[1] = (int)out[1];
  emxEnsureCapacity_creal32_T(y, i0);
  loop_ub = (int)out[0] * (int)out[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    y->data[i0].re = 0.0F;
    y->data[i0].im = 0.0F;
  }

  istart = 1.0F;
  emxInit_creal32_T(&X, 2);
  emxInit_creal32_T(&Y, 2);
  emxInit_int8_T(&r1, 1);
  emxInit_real32_T(&d_x, 2);
  emxInit_real32_T(&costab, 2);
  emxInit_real32_T(&sintab, 2);
  emxInit_real32_T(&sintabinv, 2);
  emxInit_creal32_T(&b_y, 2);
  while (istart <= b_fcnOutput) {
    y_re = istart + ((nfft - 1.0F) + 1.0F);
    y_re--;
    u1 = (float)b_fcnOutput;
    if (y_re < u1) {
      u1 = y_re;
    }

    if (u1 - istart == 0.0F) {
      i0 = r1->size[0];
      r1->size[0] = (int)nfft;
      emxEnsureCapacity_int8_T(r1, i0);
      loop_ub = (int)nfft;
      for (i0 = 0; i0 < loop_ub; i0++) {
        r1->data[i0] = 0;
      }

      loop_ub = x->size[1];
      i0 = X->size[0] * X->size[1];
      X->size[0] = r1->size[0];
      X->size[1] = loop_ub;
      emxEnsureCapacity_creal32_T(X, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_x = r1->size[0];
        for (i1 = 0; i1 < b_x; i1++) {
          X->data[i1 + X->size[0] * i0].re = x->data[((int)istart + x->size[0] *
            i0) - 1];
          X->data[i1 + X->size[0] * i0].im = 0.0F;
        }
      }

      /*  need to fft a scalar */
    } else {
      if (istart > u1) {
        i0 = 0;
        i1 = 0;
      } else {
        i0 = (int)istart - 1;
        i1 = (int)u1;
      }

      loop_ub = x->size[1];
      i2 = d_x->size[0] * d_x->size[1];
      d_x->size[0] = i1 - i0;
      d_x->size[1] = loop_ub;
      emxEnsureCapacity_real32_T(d_x, i2);
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_x = i1 - i0;
        for (i3 = 0; i3 < b_x; i3++) {
          d_x->data[i3 + d_x->size[0] * i2] = x->data[(i0 + i3) + x->size[0] *
            i2];
        }
      }

      i2 = x->size[1];
      if ((i1 - i0 == 0) || (i2 == 0)) {
        i2 = x->size[1];
        i3 = X->size[0] * X->size[1];
        X->size[0] = (int)nfft;
        X->size[1] = i2;
        emxEnsureCapacity_creal32_T(X, i3);
        if ((int)nfft > i1 - i0) {
          i0 = X->size[0] * X->size[1];
          emxEnsureCapacity_creal32_T(X, i0);
          loop_ub = X->size[1];
          for (i0 = 0; i0 < loop_ub; i0++) {
            b_x = X->size[0];
            for (i1 = 0; i1 < b_x; i1++) {
              X->data[i1 + X->size[0] * i0].re = 0.0F;
              X->data[i1 + X->size[0] * i0].im = 0.0F;
            }
          }
        }
      } else {
        useRadix2 = (((int)nfft & ((int)nfft - 1)) == 0);
        get_algo_sizes((int)nfft, useRadix2, &b_x, &unnamed_idx_1);
        generate_twiddle_tables(unnamed_idx_1, useRadix2, costab, sintab,
          sintabinv);
        if (useRadix2) {
          b_r2br_r2dit_trig(d_x, (int)nfft, costab, sintab, X);
        } else {
          dobluesteinfft(d_x, b_x, (int)nfft, costab, sintab, sintabinv, X);
        }
      }
    }

    i0 = Y->size[0] * Y->size[1];
    Y->size[0] = X->size[0];
    Y->size[1] = X->size[1];
    emxEnsureCapacity_creal32_T(Y, i0);
    loop_ub = X->size[0] * X->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      y_re = X->data[i0].re;
      u1 = X->data[i0].im;
      B_re = B->data[i0].re;
      B_im = B->data[i0].im;
      Y->data[i0].re = y_re * B_re - u1 * B_im;
      Y->data[i0].im = y_re * B_im + u1 * B_re;
    }

    b_Y[0] = Y->size[1];
    c_Y = *Y;
    c_Y.size = (int *)&b_Y;
    c_Y.numDimensions = 1;
    b_fft(&c_Y, Y->size[1], r0);
    b_x = Y->size[1];
    i0 = Y->size[0] * Y->size[1];
    Y->size[0] = 1;
    Y->size[1] = b_x;
    emxEnsureCapacity_creal32_T(Y, i0);
    for (i0 = 0; i0 < b_x; i0++) {
      Y->data[Y->size[0] * i0] = r0->data[i0];
    }

    y_re = istart + nfft;
    if ((b_fcnOutput < y_re - 1.0F) || rtIsNaNF(y_re - 1.0F)) {
      y_re = (float)b_fcnOutput;
    } else {
      y_re--;
    }

    if (istart > y_re) {
      i0 = 0;
      i1 = 0;
      i2 = 0;
    } else {
      i0 = (int)istart - 1;
      i1 = (int)y_re;
      i2 = (int)istart - 1;
    }

    b_x = y->size[1];
    i3 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = i1 - i0;
    b_y->size[1] = b_x;
    emxEnsureCapacity_creal32_T(b_y, i3);
    for (i3 = 0; i3 < b_x; i3++) {
      loop_ub = i1 - i0;
      for (unnamed_idx_1 = 0; unnamed_idx_1 < loop_ub; unnamed_idx_1++) {
        b_y->data[unnamed_idx_1 + b_y->size[0] * i3].re = y->data[(i0 +
          unnamed_idx_1) + y->size[0] * i3].re + Y->data[unnamed_idx_1 + Y->
          size[0] * i3].re;
        b_y->data[unnamed_idx_1 + b_y->size[0] * i3].im = y->data[(i0 +
          unnamed_idx_1) + y->size[0] * i3].im + Y->data[unnamed_idx_1 + Y->
          size[0] * i3].im;
      }
    }

    loop_ub = b_y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_x = b_y->size[0];
      for (i1 = 0; i1 < b_x; i1++) {
        y->data[(i2 + i1) + y->size[0] * i0] = b_y->data[i1 + b_y->size[0] * i0];
      }
    }

    istart += (nfft - 1.0F) + 1.0F;
  }

  emxFree_creal32_T(&r0);
  emxFree_real32_T(&sintabinv);
  emxFree_real32_T(&sintab);
  emxFree_real32_T(&costab);
  emxFree_real32_T(&d_x);
  emxFree_int8_T(&r1);
  emxFree_creal32_T(&Y);
  emxFree_creal32_T(&X);
  emxFree_creal32_T(&B);
  emxInit_real32_T1(&d_b, 1);
  i0 = d_b->size[0];
  d_b->size[0] = b->size[1];
  emxEnsureCapacity_real32_T1(d_b, i0);
  loop_ub = b->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_b->data[i0] = 0.0F;
  }

  if (any(d_b)) {
    useRadix2 = true;
  } else {
    i0 = d_b->size[0];
    d_b->size[0] = x->size[0] * x->size[1];
    emxEnsureCapacity_real32_T1(d_b, i0);
    loop_ub = x->size[0] * x->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      d_b->data[i0] = 0.0F;
    }

    if (any(d_b)) {
      useRadix2 = true;
    } else {
      useRadix2 = false;
    }
  }

  emxFree_real32_T(&d_b);
  if (!useRadix2) {
    loop_ub = y->size[0] * y->size[1] - 1;
    i0 = y->size[0] * y->size[1];
    emxEnsureCapacity_creal32_T(y, i0);
    for (i0 = 0; i0 <= loop_ub; i0++) {
      y_re = y->data[i0].re;
      y->data[i0].re = y_re;
      y->data[i0].im = 0.0F;
    }
  }

  if ((fcnOutput == 1.0F) && (y->size[1] == 1.0F)) {
    b_x = y->size[0];
    i0 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = b_x;
    emxEnsureCapacity_creal32_T(b_y, i0);
    for (i0 = 0; i0 < b_x; i0++) {
      b_y->data[b_y->size[0] * i0] = y->data[i0];
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = b_y->size[0];
    y->size[1] = b_y->size[1];
    emxEnsureCapacity_creal32_T(y, i0);
    loop_ub = b_y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_x = b_y->size[0];
      for (i1 = 0; i1 < b_x; i1++) {
        y->data[i1 + y->size[0] * i0] = b_y->data[i1 + b_y->size[0] * i0];
      }
    }

    /*  turn column back into a row */
  }

  emxFree_creal32_T(&b_y);
}

/* End of code generation (fftfiltering.c) */
