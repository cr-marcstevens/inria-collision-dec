#include "custom_brilliantrussian.h"
#include "m4ri/brilliantrussian.h"
#include "m4ri/misc.h"
#include "m4ri/echelonform.h"
#include "m4ri/graycode.h"


static inline int _mzd_gauss_submatrix(mzd_t *A, rci_t r, rci_t c, rci_t end_row, int k) {
  rci_t start_row = r;
  int found;
  rci_t j;
  for (j = c; j < c+k; ++j) {
    found = 0;
    for (rci_t i = start_row; i < end_row; ++i) {
      /* first we need to clear the first columns */
      for (int l = 0; l < j - c; ++l)
        if (mzd_read_bit(A, i, c+l))
          mzd_row_add_offset(A, i, r+l, c+l);

      /* pivot? */
      if (mzd_read_bit(A, i, j)) {
        mzd_row_swap(A, i, start_row);
        start_row++;
        found = 1;
        break;
      }
    }
    if (found == 0) {
      break;
    }
  }
  __M4RI_DD_MZD(A);
  __M4RI_DD_INT(j - c);
  return j - c;
}

static inline int _mzd_gauss_submatrix_top(mzd_t *A, rci_t r, rci_t c, int k) {
  rci_t start_row = r;
  for (rci_t j = c; j < c + k; ++j) {
    for (rci_t l = r; l < start_row; ++l) {
      if (mzd_read_bit(A, l, j)) {
        mzd_row_add_offset(A, l, start_row, j);
      }
    }
    ++start_row;
  }
  __M4RI_DD_MZD(A);
  __M4RI_DD_INT(k);
  return k;
}

static inline void _mzd_copy_back_rows(mzd_t *A, mzd_t const *U, rci_t r, rci_t c, int k) {
  wi_t const startblock = c / m4ri_radix;
  wi_t const width = A->width - startblock;
  for (int i = 0; i < k; ++i) {
    word const *const src = U->rows[i] + startblock;
    word *const dst = A->rows[r+i] + startblock;
    for (wi_t j = 0; j < width; ++j) {
      dst[j] = src[j];
    }
  }
  __M4RI_DD_MZD(A);
}

static inline int _mzd_gauss_submatrix_full(mzd_t *A, rci_t r, rci_t c, rci_t end_row, int k) {
  assert(k <= m4ri_radix);
  rci_t start_row = r;
  rci_t j;
  for (j = c; j < c + k; ++j) {
    int found = 0;
    for (rci_t i = start_row; i < end_row; ++i) {
      /* first we need to clear the first columns */
      word const tmp = mzd_read_bits(A, i, c, j - c + 1);
      if(tmp) {
        for (int l = 0; l < j - c; ++l)
          if (__M4RI_GET_BIT(tmp, l))
            mzd_row_add_offset(A, i, r+l, c+l);

        /* pivot? */
        if (mzd_read_bit(A, i, j)) {
          mzd_row_swap(A, i, start_row);
          /* clear above */
          for (rci_t l = r; l < start_row; ++l) {
            if (mzd_read_bit(A, l, j)) {
              mzd_row_add_offset(A, l, start_row, j);
            }
          }
          ++start_row;
          found = 1;
          break;
        }
      }
    }
    if (found == 0) {
      break;
    }
  }
  __M4RI_DD_MZD(A);
  __M4RI_DD_INT(j - c);
  return j - c;
}

rci_t _mzd_partial_echelonize_m4ri(mzd_t *A, int const full, int k, int heuristic, double const threshold, rci_t nb_col) {
  /**
	 * same as _mzd_echelonize_m4ri but stops after nb_col columns 
   */
  rci_t const ncols = A->ncols;

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.75 * __M4RI_TWOPOW(k) * ncols > __M4RI_CPU_L2_CACHE / 2.0)
      k -= 1;
  }
  int kk = 6 * k;

  mzd_t *U  = mzd_init(kk, ncols);
  mzd_t *T0 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T1 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T2 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T3 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T4 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T5 = mzd_init(__M4RI_TWOPOW(k), ncols);
  rci_t *L0 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));

  rci_t last_check = 0;
  rci_t r = 0;
  rci_t c = 0;

  if (heuristic) {
    if (c < ncols && r < A->nrows && _mzd_density(A, 32, 0, 0) >= threshold) {
      wi_t const tmp = c / m4ri_radix;
      rci_t const tmp2 = tmp * m4ri_radix;
      mzd_t *Abar = mzd_init_window(A, r, tmp2, A->nrows, ncols);
      r += mzd_echelonize_pluq(Abar, full);
      mzd_free(Abar);
      c = ncols;
    }
  }

  while(c < nb_col) {
    if (heuristic && c > (last_check + 256)) {
      last_check = c;
      if (c < ncols && r < A->nrows && _mzd_density(A, 32, r, c) >= threshold) {
        mzd_t *Abar = mzd_init_window(A, r, (c / m4ri_radix) * m4ri_radix, A->nrows, ncols);
        if (!full) {
          r += mzd_echelonize_pluq(Abar, full);
        } else {
          rci_t r2 = mzd_echelonize_pluq(Abar, full);
          if (r > 0)
            _mzd_top_echelonize_m4ri(A, 0, r, c, r);
          r += r2;
        }
        mzd_free(Abar);
        break;
      }
    }

    if(c + kk > nb_col) {
      kk = nb_col - c;
    }
    int kbar;
    if (full) {
      kbar = _mzd_gauss_submatrix_full(A, r, c, A->nrows, kk);
    } else {
      kbar = _mzd_gauss_submatrix(A, r, c, A->nrows, kk);
      /* this isn't necessary, adapt make_table */
      U = mzd_submatrix(U, A, r, 0, r + kbar, ncols);
      _mzd_gauss_submatrix_top(A, r, c, kbar);
    }

    if (kbar > 5 * k) {
      int const rem = kbar % 6;
      int const ka = kbar / 6 + ((rem >= 5) ? 1 : 0);
      int const kb = kbar / 6 + ((rem >= 4) ? 1 : 0);
      int const kc = kbar / 6 + ((rem >= 3) ? 1 : 0);
      int const kd = kbar / 6 + ((rem >= 2) ? 1 : 0);
      int const ke = kbar / 6 + ((rem >= 1) ? 1 : 0);;
      int const kf = kbar / 6;

      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
        mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
        mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
        mzd_make_table(A, r+ka+kb+kc+kd+ke, c, kf, T5, L5);
      }
      if(kbar == kk)
        mzd_process_rows6(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);
      if(full)
        mzd_process_rows6(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);

    } else if (kbar > 4 * k) {
      int const rem = kbar % 5;
      int const ka = kbar / 5 + ((rem >= 4) ? 1 : 0);
      int const kb = kbar / 5 + ((rem >= 3) ? 1 : 0);
      int const kc = kbar / 5 + ((rem >= 2) ? 1 : 0);
      int const kd = kbar / 5 + ((rem >= 1) ? 1 : 0);
      int const ke = kbar / 5;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
        mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
        mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
      }
      if(kbar == kk)
        mzd_process_rows5(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      if(full)
        mzd_process_rows5(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);

    } else if (kbar > 3 * k) {
      int const rem = kbar % 4;
      int const ka = kbar / 4 + ((rem >= 3) ? 1 : 0);
      int const kb = kbar / 4 + ((rem >= 2) ? 1 : 0);
      int const kc = kbar / 4 + ((rem >= 1) ? 1 : 0);
      int const kd = kbar / 4;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
        mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      }
      if(kbar == kk)
        mzd_process_rows4(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      if(full)
        mzd_process_rows4(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);

    } else if (kbar > 2 * k) {
      int const rem = kbar % 3;
      int const ka = kbar / 3 + ((rem >= 2) ? 1 : 0);
      int const kb = kbar / 3 + ((rem >= 1) ? 1 : 0);
      int const kc = kbar / 3;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      }
      if(kbar == kk)
        mzd_process_rows3(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2);
      if(full)
        mzd_process_rows3(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2);

    } else if (kbar > k) {
      int const ka = kbar / 2;
      int const kb = kbar - ka;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
      }
      if(kbar == kk)
        mzd_process_rows2(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1);
      if(full)
        mzd_process_rows2(A, 0, r, c, kbar, T0, L0, T1, L1);

    } else if(kbar > 0) {
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, kbar, T0, L0);
      }
      if(kbar == kk)
        mzd_process_rows(A, r+kbar, A->nrows, c, kbar, T0, L0);
      if(full)
        mzd_process_rows(A, 0, r, c, kbar, T0, L0);
    }

    if (!full) {
      _mzd_copy_back_rows(A, U, r, c, kbar);
    }

    r += kbar;
    c += kbar;
    if(kk != kbar) {
      rci_t cbar;
      rci_t rbar;
      if (mzd_find_pivot(A, r, c, &rbar, &cbar)) {
        c = cbar;
        mzd_row_swap(A, r, rbar);
      } else {
        break;
      }
      //c++;
    }
  }

  mzd_free(T0);
  m4ri_mm_free(L0);
  mzd_free(T1);
  m4ri_mm_free(L1);
  mzd_free(T2);
  m4ri_mm_free(L2);
  mzd_free(T3);
  m4ri_mm_free(L3);
  mzd_free(T4);
  m4ri_mm_free(L4);
  mzd_free(T5);
  m4ri_mm_free(L5);
  mzd_free(U);

  __M4RI_DD_MZD(A);
  __M4RI_DD_RCI(r);
  return r;
}

