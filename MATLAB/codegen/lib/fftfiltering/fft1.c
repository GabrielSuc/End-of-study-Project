/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fft1.c
 *
 * Code generation for function 'fft1'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "fftfiltering.h"
#include "fft1.h"
#include "fftfiltering_emxutil.h"
#include "bluesteinSetup.h"

/* Function Declarations */
static void b_r2br_r2dit_trig_impl(const emxArray_creal32_T *x, int
  unsigned_nRows, const emxArray_real32_T *costab, const emxArray_real32_T
  *sintab, emxArray_creal32_T *y);
static void r2br_r2dit_trig(const emxArray_creal32_T *x, int n1_unsigned, const
  emxArray_real32_T *costab, const emxArray_real32_T *sintab, emxArray_creal32_T
  *y);
static void r2br_r2dit_trig_impl(const emxArray_creal32_T *x, int unsigned_nRows,
  const emxArray_real32_T *costab, const emxArray_real32_T *sintab,
  emxArray_creal32_T *y);

/* Function Definitions */
static void b_r2br_r2dit_trig_impl(const emxArray_creal32_T *x, int
  unsigned_nRows, const emxArray_real32_T *costab, const emxArray_real32_T
  *sintab, emxArray_creal32_T *y)
{
  int j;
  int nRowsD2;
  int nRowsD4;
  int iy;
  int iDelta;
  int ix;
  int ju;
  int i;
  boolean_T tst;
  float temp_re;
  float temp_im;
  float twid_re;
  float twid_im;
  int ihi;
  j = x->size[0];
  if (!(j < unsigned_nRows)) {
    j = unsigned_nRows;
  }

  nRowsD2 = unsigned_nRows / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity_creal32_T1(y, iy);
  if (unsigned_nRows > x->size[0]) {
    iDelta = y->size[0];
    iy = y->size[0];
    y->size[0] = iDelta;
    emxEnsureCapacity_creal32_T1(y, iy);
    for (iy = 0; iy < iDelta; iy++) {
      y->data[iy].re = 0.0F;
      y->data[iy].im = 0.0F;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy] = x->data[ix];
    iDelta = unsigned_nRows;
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy] = x->data[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= unsigned_nRows - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  iDelta = 2;
  iy = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += iy) {
      temp_re = y->data[i + iDelta].re;
      temp_im = y->data[i + iDelta].im;
      y->data[i + iDelta].re = y->data[i].re - temp_re;
      y->data[i + iDelta].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    ju = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      twid_re = costab->data[j];
      twid_im = sintab->data[j];
      i = ju;
      ihi = ju + ix;
      while (i < ihi) {
        temp_re = twid_re * y->data[i + iDelta].re - twid_im * y->data[i +
          iDelta].im;
        temp_im = twid_re * y->data[i + iDelta].im + twid_im * y->data[i +
          iDelta].re;
        y->data[i + iDelta].re = y->data[i].re - temp_re;
        y->data[i + iDelta].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += iy;
      }

      ju++;
    }

    nRowsD4 /= 2;
    iDelta = iy;
    iy += iy;
    ix -= iDelta;
  }
}

static void r2br_r2dit_trig(const emxArray_creal32_T *x, int n1_unsigned, const
  emxArray_real32_T *costab, const emxArray_real32_T *sintab, emxArray_creal32_T
  *y)
{
  float r;
  int i4;
  int loop_ub;
  b_r2br_r2dit_trig_impl(x, n1_unsigned, costab, sintab, y);
  if (y->size[0] > 1) {
    r = 1.0F / (float)y->size[0];
    i4 = y->size[0];
    emxEnsureCapacity_creal32_T1(y, i4);
    loop_ub = y->size[0];
    for (i4 = 0; i4 < loop_ub; i4++) {
      y->data[i4].re *= r;
      y->data[i4].im *= r;
    }
  }
}

static void r2br_r2dit_trig_impl(const emxArray_creal32_T *x, int unsigned_nRows,
  const emxArray_real32_T *costab, const emxArray_real32_T *sintab,
  emxArray_creal32_T *y)
{
  int j;
  int nRowsD2;
  int nRowsD4;
  int iy;
  int iDelta;
  int ix;
  int ju;
  int i;
  boolean_T tst;
  float temp_re;
  float temp_im;
  float twid_re;
  float twid_im;
  int ihi;
  j = x->size[0];
  if (!(j < unsigned_nRows)) {
    j = unsigned_nRows;
  }

  nRowsD2 = unsigned_nRows / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity_creal32_T1(y, iy);
  if (unsigned_nRows > x->size[0]) {
    iDelta = y->size[0];
    iy = y->size[0];
    y->size[0] = iDelta;
    emxEnsureCapacity_creal32_T1(y, iy);
    for (iy = 0; iy < iDelta; iy++) {
      y->data[iy].re = 0.0F;
      y->data[iy].im = 0.0F;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy] = x->data[ix];
    iDelta = unsigned_nRows;
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy] = x->data[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= unsigned_nRows - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  iDelta = 2;
  iy = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += iy) {
      temp_re = y->data[i + iDelta].re;
      temp_im = y->data[i + iDelta].im;
      y->data[i + iDelta].re = y->data[i].re - temp_re;
      y->data[i + iDelta].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    ju = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      twid_re = costab->data[j];
      twid_im = sintab->data[j];
      i = ju;
      ihi = ju + ix;
      while (i < ihi) {
        temp_re = twid_re * y->data[i + iDelta].re - twid_im * y->data[i +
          iDelta].im;
        temp_im = twid_re * y->data[i + iDelta].im + twid_im * y->data[i +
          iDelta].re;
        y->data[i + iDelta].re = y->data[i].re - temp_re;
        y->data[i + iDelta].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += iy;
      }

      ju++;
    }

    nRowsD4 /= 2;
    iDelta = iy;
    iy += iy;
    ix -= iDelta;
  }
}

void b_fft(const emxArray_creal32_T *x, int n, emxArray_creal32_T *y)
{
  emxArray_real32_T *costab1q;
  emxArray_real32_T *costab;
  emxArray_real32_T *sintab;
  emxArray_real32_T *sintabinv;
  emxArray_creal32_T *wwc;
  emxArray_creal32_T *fv;
  emxArray_creal32_T *r4;
  emxArray_creal32_T *b_fv;
  boolean_T useRadix2;
  int nInt2;
  int N2blue;
  int nd2;
  float e;
  int nRowsD4;
  int k;
  int nInt2m1;
  int b_y;
  float denom_im;
  float wwc_re;
  float fv_im;
  float wwc_im;
  float fv_re;
  float b_fv_re;
  float b_fv_im;
  emxInit_real32_T(&costab1q, 2);
  emxInit_real32_T(&costab, 2);
  emxInit_real32_T(&sintab, 2);
  emxInit_real32_T(&sintabinv, 2);
  emxInit_creal32_T1(&wwc, 1);
  emxInit_creal32_T1(&fv, 1);
  emxInit_creal32_T1(&r4, 1);
  emxInit_creal32_T1(&b_fv, 1);
  if ((x->size[0] == 0) || (n == 0)) {
    nInt2 = y->size[0];
    y->size[0] = n;
    emxEnsureCapacity_creal32_T1(y, nInt2);
    if (n > x->size[0]) {
      nd2 = y->size[0];
      nInt2 = y->size[0];
      y->size[0] = nd2;
      emxEnsureCapacity_creal32_T1(y, nInt2);
      for (nInt2 = 0; nInt2 < nd2; nInt2++) {
        y->data[nInt2].re = 0.0F;
        y->data[nInt2].im = 0.0F;
      }
    }
  } else {
    useRadix2 = ((n & (n - 1)) == 0);
    get_algo_sizes(n, useRadix2, &N2blue, &nd2);
    e = 6.28318548F / (float)nd2;
    nRowsD4 = nd2 / 2 / 2;
    nInt2 = costab1q->size[0] * costab1q->size[1];
    costab1q->size[0] = 1;
    costab1q->size[1] = nRowsD4 + 1;
    emxEnsureCapacity_real32_T(costab1q, nInt2);
    costab1q->data[0] = 1.0F;
    nd2 = nRowsD4 / 2;
    for (k = 1; k <= nd2; k++) {
      costab1q->data[k] = (float)cos(e * (float)k);
    }

    for (k = nd2 + 1; k < nRowsD4; k++) {
      costab1q->data[k] = (float)sin(e * (float)(nRowsD4 - k));
    }

    costab1q->data[nRowsD4] = 0.0F;
    if (!useRadix2) {
      nd2 = costab1q->size[1] - 1;
      nRowsD4 = (costab1q->size[1] - 1) << 1;
      nInt2 = costab->size[0] * costab->size[1];
      costab->size[0] = 1;
      costab->size[1] = nRowsD4 + 1;
      emxEnsureCapacity_real32_T(costab, nInt2);
      nInt2 = sintab->size[0] * sintab->size[1];
      sintab->size[0] = 1;
      sintab->size[1] = nRowsD4 + 1;
      emxEnsureCapacity_real32_T(sintab, nInt2);
      costab->data[0] = 1.0F;
      sintab->data[0] = 0.0F;
      nInt2 = sintabinv->size[0] * sintabinv->size[1];
      sintabinv->size[0] = 1;
      sintabinv->size[1] = nRowsD4 + 1;
      emxEnsureCapacity_real32_T(sintabinv, nInt2);
      for (k = 1; k <= nd2; k++) {
        sintabinv->data[k] = costab1q->data[nd2 - k];
      }

      for (k = costab1q->size[1]; k <= nRowsD4; k++) {
        sintabinv->data[k] = costab1q->data[k - nd2];
      }

      for (k = 1; k <= nd2; k++) {
        costab->data[k] = costab1q->data[k];
        sintab->data[k] = -costab1q->data[nd2 - k];
      }

      for (k = costab1q->size[1]; k <= nRowsD4; k++) {
        costab->data[k] = -costab1q->data[nRowsD4 - k];
        sintab->data[k] = -costab1q->data[k - nd2];
      }
    } else {
      nd2 = costab1q->size[1] - 1;
      nRowsD4 = (costab1q->size[1] - 1) << 1;
      nInt2 = costab->size[0] * costab->size[1];
      costab->size[0] = 1;
      costab->size[1] = nRowsD4 + 1;
      emxEnsureCapacity_real32_T(costab, nInt2);
      nInt2 = sintab->size[0] * sintab->size[1];
      sintab->size[0] = 1;
      sintab->size[1] = nRowsD4 + 1;
      emxEnsureCapacity_real32_T(sintab, nInt2);
      costab->data[0] = 1.0F;
      sintab->data[0] = 0.0F;
      for (k = 1; k <= nd2; k++) {
        costab->data[k] = costab1q->data[k];
        sintab->data[k] = costab1q->data[nd2 - k];
      }

      for (k = costab1q->size[1]; k <= nRowsD4; k++) {
        costab->data[k] = -costab1q->data[nRowsD4 - k];
        sintab->data[k] = costab1q->data[k - nd2];
      }

      nInt2 = sintabinv->size[0] * sintabinv->size[1];
      sintabinv->size[0] = 1;
      sintabinv->size[1] = 0;
      emxEnsureCapacity_real32_T(sintabinv, nInt2);
    }

    if (useRadix2) {
      r2br_r2dit_trig(x, n, costab, sintab, y);
    } else {
      nInt2m1 = (n + n) - 1;
      nInt2 = wwc->size[0];
      wwc->size[0] = nInt2m1;
      emxEnsureCapacity_creal32_T1(wwc, nInt2);
      nd2 = n;
      nRowsD4 = 0;
      wwc->data[n - 1].re = 1.0F;
      wwc->data[n - 1].im = 0.0F;
      nInt2 = n << 1;
      for (k = 1; k < n; k++) {
        b_y = (k << 1) - 1;
        if (nInt2 - nRowsD4 <= b_y) {
          nRowsD4 += b_y - nInt2;
        } else {
          nRowsD4 += b_y;
        }

        denom_im = 3.14159274F * (float)nRowsD4 / (float)n;
        if (denom_im == 0.0F) {
          e = 1.0F;
          denom_im = 0.0F;
        } else {
          e = (float)cos(denom_im);
          denom_im = (float)sin(denom_im);
        }

        wwc->data[nd2 - 2].re = e;
        wwc->data[nd2 - 2].im = -denom_im;
        nd2--;
      }

      nd2 = 0;
      for (k = nInt2m1 - 1; k >= n; k--) {
        wwc->data[k] = wwc->data[nd2];
        nd2++;
      }

      nRowsD4 = x->size[0];
      if (n < nRowsD4) {
        nRowsD4 = n;
      }

      nInt2 = y->size[0];
      y->size[0] = n;
      emxEnsureCapacity_creal32_T1(y, nInt2);
      if (n > x->size[0]) {
        nd2 = y->size[0];
        nInt2 = y->size[0];
        y->size[0] = nd2;
        emxEnsureCapacity_creal32_T1(y, nInt2);
        for (nInt2 = 0; nInt2 < nd2; nInt2++) {
          y->data[nInt2].re = 0.0F;
          y->data[nInt2].im = 0.0F;
        }
      }

      nd2 = 0;
      for (k = 0; k < nRowsD4; k++) {
        e = wwc->data[(n + k) - 1].re;
        denom_im = wwc->data[(n + k) - 1].im;
        wwc_re = x->data[nd2].re;
        fv_im = x->data[nd2].im;
        wwc_im = x->data[nd2].im;
        fv_re = x->data[nd2].re;
        y->data[k].re = e * wwc_re + denom_im * fv_im;
        y->data[k].im = e * wwc_im - denom_im * fv_re;
        nd2++;
      }

      while (nRowsD4 + 1 <= n) {
        y->data[nRowsD4].re = 0.0F;
        y->data[nRowsD4].im = 0.0F;
        nRowsD4++;
      }

      r2br_r2dit_trig_impl(y, N2blue, costab, sintab, fv);
      b_r2br_r2dit_trig_impl(wwc, N2blue, costab, sintab, r4);
      nInt2 = b_fv->size[0];
      b_fv->size[0] = fv->size[0];
      emxEnsureCapacity_creal32_T1(b_fv, nInt2);
      nd2 = fv->size[0];
      for (nInt2 = 0; nInt2 < nd2; nInt2++) {
        b_fv_re = fv->data[nInt2].re;
        b_fv_im = fv->data[nInt2].im;
        e = r4->data[nInt2].re;
        denom_im = r4->data[nInt2].im;
        b_fv->data[nInt2].re = b_fv_re * e - b_fv_im * denom_im;
        b_fv->data[nInt2].im = b_fv_re * denom_im + b_fv_im * e;
      }

      r2br_r2dit_trig(b_fv, N2blue, costab, sintabinv, fv);
      nd2 = 0;
      for (k = (int)(float)n - 1; k < wwc->size[0]; k++) {
        e = wwc->data[k].re;
        b_fv_re = fv->data[k].re;
        denom_im = wwc->data[k].im;
        b_fv_im = fv->data[k].im;
        wwc_re = wwc->data[k].re;
        fv_im = fv->data[k].im;
        wwc_im = wwc->data[k].im;
        fv_re = fv->data[k].re;
        y->data[nd2].re = e * b_fv_re + denom_im * b_fv_im;
        y->data[nd2].im = wwc_re * fv_im - wwc_im * fv_re;
        e = wwc->data[k].re;
        b_fv_re = fv->data[k].re;
        denom_im = wwc->data[k].im;
        b_fv_im = fv->data[k].im;
        wwc_re = wwc->data[k].re;
        fv_im = fv->data[k].im;
        wwc_im = wwc->data[k].im;
        fv_re = fv->data[k].re;
        y->data[nd2].re = e * b_fv_re + denom_im * b_fv_im;
        y->data[nd2].im = wwc_re * fv_im - wwc_im * fv_re;
        e = y->data[nd2].re;
        denom_im = y->data[nd2].im;
        if (denom_im == 0.0F) {
          y->data[nd2].re = e / (float)n;
          y->data[nd2].im = 0.0F;
        } else if (e == 0.0F) {
          y->data[nd2].re = 0.0F;
          y->data[nd2].im = denom_im / (float)n;
        } else {
          y->data[nd2].re = e / (float)n;
          y->data[nd2].im = denom_im / (float)n;
        }

        nd2++;
      }
    }
  }

  emxFree_creal32_T(&b_fv);
  emxFree_creal32_T(&r4);
  emxFree_creal32_T(&fv);
  emxFree_creal32_T(&wwc);
  emxFree_real32_T(&sintabinv);
  emxFree_real32_T(&sintab);
  emxFree_real32_T(&costab);
  emxFree_real32_T(&costab1q);
}

void b_r2br_r2dit_trig(const emxArray_real32_T *x, int n1_unsigned, const
  emxArray_real32_T *costab, const emxArray_real32_T *sintab, emxArray_creal32_T
  *y)
{
  emxArray_creal32_T *rwork;
  int ub_loop;
  int n1;
  unsigned int sx[2];
  int nrows;
  int sz[2];
  int k;
  int loop_ub;
  int b_loop_ub;
  int i5;
  int xoff;
  int iy;
  int j;
  int nRowsD2;
  int nRowsD4;
  int iDelta;
  int ju;
  int i;
  boolean_T tst;
  float temp_re;
  float temp_im;
  float twid_re;
  float twid_im;
  int ihi;
  emxInit_creal32_T1(&rwork, 1);
  for (ub_loop = 0; ub_loop < 2; ub_loop++) {
    sx[ub_loop] = (unsigned int)x->size[ub_loop];
  }

  n1 = n1_unsigned;
  nrows = x->size[0];
  for (ub_loop = 0; ub_loop < 2; ub_loop++) {
    sz[ub_loop] = x->size[ub_loop];
  }

  ub_loop = y->size[0] * y->size[1];
  y->size[0] = n1_unsigned;
  y->size[1] = sz[1];
  emxEnsureCapacity_creal32_T(y, ub_loop);
  if (n1_unsigned > x->size[0]) {
    ub_loop = y->size[0] * y->size[1];
    emxEnsureCapacity_creal32_T(y, ub_loop);
    loop_ub = y->size[1];
    for (ub_loop = 0; ub_loop < loop_ub; ub_loop++) {
      b_loop_ub = y->size[0];
      for (i5 = 0; i5 < b_loop_ub; i5++) {
        y->data[i5 + y->size[0] * ub_loop].re = 0.0F;
        y->data[i5 + y->size[0] * ub_loop].im = 0.0F;
      }
    }
  }

  ub_loop = (int)sx[1];

#pragma omp parallel \
 num_threads(omp_get_max_threads()) \
 private(rwork,xoff,iy,j,nRowsD2,nRowsD4,iDelta,ju,i,tst,temp_re,temp_im,twid_re,twid_im,ihi)

  {
    emxInit_creal32_T1(&rwork, 1);

#pragma omp for nowait

    for (k = 1; k <= ub_loop; k++) {
      xoff = (k - 1) * nrows;
      iy = x->size[0];
      j = n1_unsigned;
      if (iy < n1_unsigned) {
        j = iy;
      }

      nRowsD2 = n1_unsigned / 2;
      nRowsD4 = nRowsD2 / 2;
      iy = rwork->size[0];
      rwork->size[0] = n1_unsigned;
      emxEnsureCapacity_creal32_T1(rwork, iy);
      if (n1_unsigned > x->size[0]) {
        iDelta = rwork->size[0];
        iy = rwork->size[0];
        rwork->size[0] = iDelta;
        emxEnsureCapacity_creal32_T1(rwork, iy);
        for (iy = 0; iy < iDelta; iy++) {
          rwork->data[iy].re = 0.0F;
          rwork->data[iy].im = 0.0F;
        }
      }

      ju = 0;
      iy = 0;
      for (i = 1; i < j; i++) {
        rwork->data[iy].re = x->data[xoff];
        rwork->data[iy].im = 0.0F;
        iDelta = n1_unsigned;
        tst = true;
        while (tst) {
          iDelta >>= 1;
          ju ^= iDelta;
          tst = ((ju & iDelta) == 0);
        }

        iy = ju;
        xoff++;
      }

      rwork->data[iy].re = x->data[xoff];
      rwork->data[iy].im = 0.0F;
      if (n1_unsigned > 1) {
        for (i = 0; i <= n1_unsigned - 2; i += 2) {
          temp_re = rwork->data[i + 1].re;
          temp_im = rwork->data[i + 1].im;
          rwork->data[i + 1].re = rwork->data[i].re - rwork->data[i + 1].re;
          rwork->data[i + 1].im = rwork->data[i].im - rwork->data[i + 1].im;
          rwork->data[i].re += temp_re;
          rwork->data[i].im += temp_im;
        }
      }

      iDelta = 2;
      iy = 4;
      ju = 1 + ((nRowsD4 - 1) << 2);
      while (nRowsD4 > 0) {
        for (i = 0; i < ju; i += iy) {
          temp_re = rwork->data[i + iDelta].re;
          temp_im = rwork->data[i + iDelta].im;
          rwork->data[i + iDelta].re = rwork->data[i].re - temp_re;
          rwork->data[i + iDelta].im = rwork->data[i].im - temp_im;
          rwork->data[i].re += temp_re;
          rwork->data[i].im += temp_im;
        }

        xoff = 1;
        for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
          twid_re = costab->data[j];
          twid_im = sintab->data[j];
          i = xoff;
          ihi = xoff + ju;
          while (i < ihi) {
            temp_re = twid_re * rwork->data[i + iDelta].re - twid_im *
              rwork->data[i + iDelta].im;
            temp_im = twid_re * rwork->data[i + iDelta].im + twid_im *
              rwork->data[i + iDelta].re;
            rwork->data[i + iDelta].re = rwork->data[i].re - temp_re;
            rwork->data[i + iDelta].im = rwork->data[i].im - temp_im;
            rwork->data[i].re += temp_re;
            rwork->data[i].im += temp_im;
            i += iy;
          }

          xoff++;
        }

        nRowsD4 /= 2;
        iDelta = iy;
        iy += iy;
        ju -= iDelta;
      }

      for (iy = 0; iy < n1; iy++) {
        y->data[iy + y->size[0] * (k - 1)] = rwork->data[iy];
      }
    }

    emxFree_creal32_T(&rwork);
  }

  emxFree_creal32_T(&rwork);
}

void dobluesteinfft(const emxArray_real32_T *x, int N2, int n1, const
                    emxArray_real32_T *costab, const emxArray_real32_T *sintab,
                    const emxArray_real32_T *sintabinv, emxArray_creal32_T *y)
{
  emxArray_creal32_T *rwork;
  emxArray_creal32_T *fv;
  emxArray_creal32_T *r3;
  emxArray_creal32_T *b_fv;
  int ub_loop;
  emxArray_creal32_T *wwc;
  unsigned int sx[2];
  int nrows;
  int sz[2];
  int k;
  int loop_ub;
  int b_loop_ub;
  int i6;
  int xoff;
  int minNrowsNx;
  int idx;
  int c_loop_ub;
  float re;
  float im;
  float fv_re;
  float fv_im;
  float wwc_re;
  float b_fv_im;
  float wwc_im;
  float b_fv_re;
  emxInit_creal32_T1(&rwork, 1);
  emxInit_creal32_T1(&fv, 1);
  emxInit_creal32_T1(&r3, 1);
  emxInit_creal32_T1(&b_fv, 1);
  for (ub_loop = 0; ub_loop < 2; ub_loop++) {
    sx[ub_loop] = (unsigned int)x->size[ub_loop];
  }

  emxInit_creal32_T1(&wwc, 1);
  bluesteinSetup(n1, wwc);
  nrows = x->size[0];
  for (ub_loop = 0; ub_loop < 2; ub_loop++) {
    sz[ub_loop] = x->size[ub_loop];
  }

  ub_loop = y->size[0] * y->size[1];
  y->size[0] = n1;
  y->size[1] = sz[1];
  emxEnsureCapacity_creal32_T(y, ub_loop);
  if (n1 > x->size[0]) {
    ub_loop = y->size[0] * y->size[1];
    emxEnsureCapacity_creal32_T(y, ub_loop);
    loop_ub = y->size[1];
    for (ub_loop = 0; ub_loop < loop_ub; ub_loop++) {
      b_loop_ub = y->size[0];
      for (i6 = 0; i6 < b_loop_ub; i6++) {
        y->data[i6 + y->size[0] * ub_loop].re = 0.0F;
        y->data[i6 + y->size[0] * ub_loop].im = 0.0F;
      }
    }
  }

  ub_loop = (int)sx[1];

#pragma omp parallel \
 num_threads(omp_get_max_threads()) \
 private(b_fv,r3,fv,rwork,xoff,minNrowsNx,idx,c_loop_ub,re,im,fv_re,fv_im,wwc_re,b_fv_im,wwc_im,b_fv_re)

  {
    emxInit_creal32_T1(&b_fv, 1);
    emxInit_creal32_T1(&r3, 1);
    emxInit_creal32_T1(&fv, 1);
    emxInit_creal32_T1(&rwork, 1);

#pragma omp for nowait

    for (k = 1; k <= ub_loop; k++) {
      xoff = (k - 1) * nrows;
      minNrowsNx = x->size[0];
      if (n1 < minNrowsNx) {
        minNrowsNx = n1;
      }

      idx = rwork->size[0];
      rwork->size[0] = n1;
      emxEnsureCapacity_creal32_T1(rwork, idx);
      if (n1 > x->size[0]) {
        c_loop_ub = rwork->size[0];
        idx = rwork->size[0];
        rwork->size[0] = c_loop_ub;
        emxEnsureCapacity_creal32_T1(rwork, idx);
        for (idx = 0; idx < c_loop_ub; idx++) {
          rwork->data[idx].re = 0.0F;
          rwork->data[idx].im = 0.0F;
        }
      }

      for (c_loop_ub = 0; c_loop_ub < minNrowsNx; c_loop_ub++) {
        re = wwc->data[(n1 + c_loop_ub) - 1].re;
        im = wwc->data[(n1 + c_loop_ub) - 1].im;
        rwork->data[c_loop_ub].re = re * x->data[xoff];
        rwork->data[c_loop_ub].im = im * -x->data[xoff];
        xoff++;
      }

      while (minNrowsNx + 1 <= n1) {
        rwork->data[minNrowsNx].re = 0.0F;
        rwork->data[minNrowsNx].im = 0.0F;
        minNrowsNx++;
      }

      r2br_r2dit_trig_impl(rwork, N2, costab, sintab, fv);
      b_r2br_r2dit_trig_impl(wwc, N2, costab, sintab, r3);
      idx = b_fv->size[0];
      b_fv->size[0] = fv->size[0];
      emxEnsureCapacity_creal32_T1(b_fv, idx);
      c_loop_ub = fv->size[0];
      for (idx = 0; idx < c_loop_ub; idx++) {
        fv_re = fv->data[idx].re;
        fv_im = fv->data[idx].im;
        re = r3->data[idx].re;
        im = r3->data[idx].im;
        b_fv->data[idx].re = fv_re * re - fv_im * im;
        b_fv->data[idx].im = fv_re * im + fv_im * re;
      }

      r2br_r2dit_trig(b_fv, N2, costab, sintabinv, fv);
      idx = 0;
      for (c_loop_ub = (int)(float)n1 - 1; c_loop_ub < wwc->size[0]; c_loop_ub++)
      {
        re = wwc->data[c_loop_ub].re;
        fv_re = fv->data[c_loop_ub].re;
        im = wwc->data[c_loop_ub].im;
        fv_im = fv->data[c_loop_ub].im;
        wwc_re = wwc->data[c_loop_ub].re;
        b_fv_im = fv->data[c_loop_ub].im;
        wwc_im = wwc->data[c_loop_ub].im;
        b_fv_re = fv->data[c_loop_ub].re;
        rwork->data[idx].re = re * fv_re + im * fv_im;
        rwork->data[idx].im = wwc_re * b_fv_im - wwc_im * b_fv_re;
        idx++;
      }

      for (idx = 0; idx < n1; idx++) {
        y->data[idx + y->size[0] * (k - 1)] = rwork->data[idx];
      }
    }

    emxFree_creal32_T(&rwork);
    emxFree_creal32_T(&fv);
    emxFree_creal32_T(&r3);
    emxFree_creal32_T(&b_fv);
  }

  emxFree_creal32_T(&wwc);
  emxFree_creal32_T(&b_fv);
  emxFree_creal32_T(&r3);
  emxFree_creal32_T(&fv);
  emxFree_creal32_T(&rwork);
}

void fft(const emxArray_real32_T *x, int n, emxArray_creal32_T *y)
{
  emxArray_real32_T *costab;
  int xidx;
  emxArray_real32_T *sintab;
  emxArray_real32_T *sintabinv;
  boolean_T useRadix2;
  int ju;
  int iheight;
  emxArray_creal32_T *wwc;
  int j;
  int minNrowsNx;
  int nRowsD2;
  int nRowsD4;
  int i;
  float twid_re;
  float twid_im;
  emxArray_creal32_T *fv;
  emxArray_creal32_T *r2;
  emxArray_creal32_T *b_fv;
  float temp_re;
  float temp_im;
  float fv_re;
  float fv_im;
  int ihi;
  float wwc_im;
  float b_fv_re;
  if (x->size[0] == 0) {
    xidx = y->size[0];
    y->size[0] = n;
    emxEnsureCapacity_creal32_T1(y, xidx);
    if (n > 0) {
      iheight = y->size[0];
      xidx = y->size[0];
      y->size[0] = iheight;
      emxEnsureCapacity_creal32_T1(y, xidx);
      for (xidx = 0; xidx < iheight; xidx++) {
        y->data[xidx].re = 0.0F;
        y->data[xidx].im = 0.0F;
      }
    }
  } else {
    emxInit_real32_T(&costab, 2);
    emxInit_real32_T(&sintab, 2);
    emxInit_real32_T(&sintabinv, 2);
    useRadix2 = ((n & (n - 1)) == 0);
    get_algo_sizes(n, useRadix2, &ju, &xidx);
    generate_twiddle_tables(xidx, useRadix2, costab, sintab, sintabinv);
    if (useRadix2) {
      j = x->size[0];
      if (!(j < n)) {
        j = n;
      }

      nRowsD2 = n / 2;
      nRowsD4 = nRowsD2 / 2;
      xidx = y->size[0];
      y->size[0] = n;
      emxEnsureCapacity_creal32_T1(y, xidx);
      if (n > x->size[0]) {
        iheight = y->size[0];
        xidx = y->size[0];
        y->size[0] = iheight;
        emxEnsureCapacity_creal32_T1(y, xidx);
        for (xidx = 0; xidx < iheight; xidx++) {
          y->data[xidx].re = 0.0F;
          y->data[xidx].im = 0.0F;
        }
      }

      minNrowsNx = 0;
      ju = 0;
      xidx = 0;
      for (i = 1; i < j; i++) {
        y->data[xidx].re = x->data[minNrowsNx];
        y->data[xidx].im = 0.0F;
        xidx = n;
        useRadix2 = true;
        while (useRadix2) {
          xidx >>= 1;
          ju ^= xidx;
          useRadix2 = ((ju & xidx) == 0);
        }

        xidx = ju;
        minNrowsNx++;
      }

      y->data[xidx].re = x->data[minNrowsNx];
      y->data[xidx].im = 0.0F;
      if (n > 1) {
        for (i = 0; i <= n - 2; i += 2) {
          temp_re = y->data[i + 1].re;
          temp_im = y->data[i + 1].im;
          y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
          y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
        }
      }

      xidx = 2;
      minNrowsNx = 4;
      iheight = 1 + ((nRowsD4 - 1) << 2);
      while (nRowsD4 > 0) {
        for (i = 0; i < iheight; i += minNrowsNx) {
          temp_re = y->data[i + xidx].re;
          temp_im = y->data[i + xidx].im;
          y->data[i + xidx].re = y->data[i].re - temp_re;
          y->data[i + xidx].im = y->data[i].im - temp_im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
        }

        ju = 1;
        for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
          twid_re = costab->data[j];
          twid_im = sintab->data[j];
          i = ju;
          ihi = ju + iheight;
          while (i < ihi) {
            temp_re = twid_re * y->data[i + xidx].re - twid_im * y->data[i +
              xidx].im;
            temp_im = twid_re * y->data[i + xidx].im + twid_im * y->data[i +
              xidx].re;
            y->data[i + xidx].re = y->data[i].re - temp_re;
            y->data[i + xidx].im = y->data[i].im - temp_im;
            y->data[i].re += temp_re;
            y->data[i].im += temp_im;
            i += minNrowsNx;
          }

          ju++;
        }

        nRowsD4 /= 2;
        xidx = minNrowsNx;
        minNrowsNx += minNrowsNx;
        iheight -= xidx;
      }
    } else {
      emxInit_creal32_T1(&wwc, 1);
      bluesteinSetup(n, wwc);
      minNrowsNx = x->size[0];
      if (n < minNrowsNx) {
        minNrowsNx = n;
      }

      xidx = y->size[0];
      y->size[0] = n;
      emxEnsureCapacity_creal32_T1(y, xidx);
      if (n > x->size[0]) {
        iheight = y->size[0];
        xidx = y->size[0];
        y->size[0] = iheight;
        emxEnsureCapacity_creal32_T1(y, xidx);
        for (xidx = 0; xidx < iheight; xidx++) {
          y->data[xidx].re = 0.0F;
          y->data[xidx].im = 0.0F;
        }
      }

      xidx = 0;
      for (iheight = 0; iheight < minNrowsNx; iheight++) {
        twid_re = wwc->data[(n + iheight) - 1].re;
        twid_im = wwc->data[(n + iheight) - 1].im;
        y->data[iheight].re = twid_re * x->data[xidx];
        y->data[iheight].im = twid_im * -x->data[xidx];
        xidx++;
      }

      while (minNrowsNx + 1 <= n) {
        y->data[minNrowsNx].re = 0.0F;
        y->data[minNrowsNx].im = 0.0F;
        minNrowsNx++;
      }

      emxInit_creal32_T1(&fv, 1);
      emxInit_creal32_T1(&r2, 1);
      emxInit_creal32_T1(&b_fv, 1);
      r2br_r2dit_trig_impl(y, ju, costab, sintab, fv);
      b_r2br_r2dit_trig_impl(wwc, ju, costab, sintab, r2);
      xidx = b_fv->size[0];
      b_fv->size[0] = fv->size[0];
      emxEnsureCapacity_creal32_T1(b_fv, xidx);
      iheight = fv->size[0];
      for (xidx = 0; xidx < iheight; xidx++) {
        fv_re = fv->data[xidx].re;
        fv_im = fv->data[xidx].im;
        twid_re = r2->data[xidx].re;
        twid_im = r2->data[xidx].im;
        b_fv->data[xidx].re = fv_re * twid_re - fv_im * twid_im;
        b_fv->data[xidx].im = fv_re * twid_im + fv_im * twid_re;
      }

      emxFree_creal32_T(&r2);
      r2br_r2dit_trig(b_fv, ju, costab, sintabinv, fv);
      xidx = 0;
      iheight = (int)(float)n - 1;
      emxFree_creal32_T(&b_fv);
      while (iheight + 1 <= wwc->size[0]) {
        twid_re = wwc->data[iheight].re;
        fv_re = fv->data[iheight].re;
        twid_im = wwc->data[iheight].im;
        fv_im = fv->data[iheight].im;
        temp_re = wwc->data[iheight].re;
        temp_im = fv->data[iheight].im;
        wwc_im = wwc->data[iheight].im;
        b_fv_re = fv->data[iheight].re;
        y->data[xidx].re = twid_re * fv_re + twid_im * fv_im;
        y->data[xidx].im = temp_re * temp_im - wwc_im * b_fv_re;
        xidx++;
        iheight++;
      }

      emxFree_creal32_T(&fv);
      emxFree_creal32_T(&wwc);
    }

    emxFree_real32_T(&sintabinv);
    emxFree_real32_T(&sintab);
    emxFree_real32_T(&costab);
  }
}

void generate_twiddle_tables(int nRows, boolean_T useRadix2, emxArray_real32_T
  *costab, emxArray_real32_T *sintab, emxArray_real32_T *sintabinv)
{
  emxArray_real32_T *costab1q;
  float e;
  int nRowsD4;
  int nd2;
  int k;
  int n2;
  emxInit_real32_T(&costab1q, 2);
  e = 6.28318548F / (float)nRows;
  nRowsD4 = nRows / 2 / 2;
  nd2 = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = nRowsD4 + 1;
  emxEnsureCapacity_real32_T(costab1q, nd2);
  costab1q->data[0] = 1.0F;
  nd2 = nRowsD4 / 2;
  for (k = 1; k <= nd2; k++) {
    costab1q->data[k] = (float)cos(e * (float)k);
  }

  for (k = nd2 + 1; k < nRowsD4; k++) {
    costab1q->data[k] = (float)sin(e * (float)(nRowsD4 - k));
  }

  costab1q->data[nRowsD4] = 0.0F;
  if (!useRadix2) {
    nRowsD4 = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity_real32_T(costab, nd2);
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity_real32_T(sintab, nd2);
    costab->data[0] = 1.0F;
    sintab->data[0] = 0.0F;
    nd2 = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = n2 + 1;
    emxEnsureCapacity_real32_T(sintabinv, nd2);
    for (k = 1; k <= nRowsD4; k++) {
      sintabinv->data[k] = costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      sintabinv->data[k] = costab1q->data[k - nRowsD4];
    }

    for (k = 1; k <= nRowsD4; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - nRowsD4];
    }
  } else {
    nRowsD4 = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity_real32_T(costab, nd2);
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity_real32_T(sintab, nd2);
    costab->data[0] = 1.0F;
    sintab->data[0] = 0.0F;
    for (k = 1; k <= nRowsD4; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - nRowsD4];
    }

    nd2 = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = 0;
    emxEnsureCapacity_real32_T(sintabinv, nd2);
  }

  emxFree_real32_T(&costab1q);
}

void get_algo_sizes(int n1, boolean_T useRadix2, int *N2blue, int *nRows)
{
  int nn1m1;
  int pmax;
  int pmin;
  boolean_T exitg1;
  int p;
  int pow2p;
  *N2blue = 1;
  if (useRadix2) {
    *nRows = n1;
  } else {
    nn1m1 = (n1 + n1) - 1;
    pmax = 31;
    if (nn1m1 <= 1) {
      pmax = 0;
    } else {
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        p = (pmin + pmax) >> 1;
        pow2p = 1 << p;
        if (pow2p == nn1m1) {
          pmax = p;
          exitg1 = true;
        } else if (pow2p > nn1m1) {
          pmax = p;
        } else {
          pmin = p;
        }
      }
    }

    *N2blue = 1 << pmax;
    *nRows = *N2blue;
  }
}

/* End of code generation (fft1.c) */
