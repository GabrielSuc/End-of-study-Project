/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bluesteinSetup.c
 *
 * Code generation for function 'bluesteinSetup'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "fftfiltering.h"
#include "bluesteinSetup.h"
#include "fftfiltering_emxutil.h"

/* Function Definitions */
void bluesteinSetup(int nRows, emxArray_creal32_T *wwc)
{
  int nInt2m1;
  int idx;
  int rt;
  int nInt2;
  int k;
  int y;
  float nt_im;
  float nt_re;
  nInt2m1 = (nRows + nRows) - 1;
  idx = wwc->size[0];
  wwc->size[0] = nInt2m1;
  emxEnsureCapacity_creal32_T1(wwc, idx);
  idx = nRows;
  rt = 0;
  wwc->data[nRows - 1].re = 1.0F;
  wwc->data[nRows - 1].im = 0.0F;
  nInt2 = nRows << 1;
  for (k = 1; k < nRows; k++) {
    y = (k << 1) - 1;
    if (nInt2 - rt <= y) {
      rt += y - nInt2;
    } else {
      rt += y;
    }

    nt_im = -3.14159274F * (float)rt / (float)nRows;
    if (nt_im == 0.0F) {
      nt_re = 1.0F;
      nt_im = 0.0F;
    } else {
      nt_re = (float)cos(nt_im);
      nt_im = (float)sin(nt_im);
    }

    wwc->data[idx - 2].re = nt_re;
    wwc->data[idx - 2].im = -nt_im;
    idx--;
  }

  idx = 0;
  for (k = nInt2m1 - 1; k >= nRows; k--) {
    wwc->data[k] = wwc->data[idx];
    idx++;
  }
}

/* End of code generation (bluesteinSetup.c) */
