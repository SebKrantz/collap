/*
 This code is adapted from the kit package: https://github.com/2005m/kit
 and licensed under a GPL-3.0 license.
*/

#include "kit.h"

#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT

// ****************************************
// This function groups a single vector
// ****************************************
SEXP dupVecIndex(SEXP x) {
  const int n = length(x);
  int K = 0, tx = TYPEOF(x), x_min = INT_MAX, x_max = INT_MIN, anyNA = 0;
  size_t M;
  // if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
    bigint:;
    const size_t n2 = 2U * (size_t) n;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  } else if(tx == INTSXP) {
    if(isFactor(x) || inherits(x, "qG")) {
      tx = 1000;
      M = isFactor(x) ? (size_t)nlevels(x) + 2 : (size_t)asInteger(getAttrib(x, install("N.groups"))) + 2;
      anyNA = !inherits(x, "na.included");
    } else {
      int *restrict p = INTEGER(x);
      // Old:
      if(n < 10 || NOGE(p[0], n) || NOGE(p[n/2], n) || NOGE(p[n-1], n)) {
        // This loop is highly optimized...
        for(int i = 0, x_tmp; i != n; ++i) {
          x_tmp = p[i];
          if(x_tmp > x_max) x_max = x_tmp;
          if(x_tmp < x_min) {
            if(x_tmp == NA_INTEGER) anyNA = 1;
            else x_min = x_tmp;
          }
        }
        double x_diff = (double)x_max - x_min;
        if(x_diff >= INT_MAX || x_diff <= INT_MIN) goto bigint; // To avoid overflows (UBSAN errors)
        x_max -= x_min;
        if(++x_max > 3 * n) goto bigint;
        M = (size_t)(x_max + 2);
        if(x_min == 0 || x_min == 1) tx = 1000;
        else x_max = NA_INTEGER;
      } else M = (size_t)n;
    }
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  SEXP ans_i = PROTECT(allocVector(INTSXP, n));
  int *restrict pans_i = INTEGER(ans_i), g = 0;
  size_t id = 0;
  switch (tx) {
  case LGLSXP:
  case 1000: // This is for factors or logical vectors where the size of the table is known
  {
    const int *restrict px = INTEGER(x);
    if(tx == 1000 && !anyNA) {
      for(int i = 0, j; i != n; ++i) {
        j = px[i];
        if(h[j]) pans_i[i] = h[j];
        else pans_i[i] = h[j] = ++g;
      }
    } else {
      for(int i = 0, j, k = (int)M-1; i != n; ++i) {
        j = (px[i] == NA_INTEGER) ? k : px[i];
        if(h[j]) pans_i[i] = h[j];
        else pans_i[i] = h[j] = ++g;
      }
    }
  } break;
  case INTSXP: {
    const int *restrict px = INTEGER(x);
    // Old:
    if(x_max == INT_MIN && M == (size_t)n) { // Faster version based on division hash...
      unsigned int iid = 0, nu = (unsigned)n;
      for (int i = 0; i != n; ++i) {
        iid = (unsigned)px[i];
        if(iid >= nu) iid %= nu;
        while(h[iid]) {
          if(px[h[iid]-1] == px[i]) {
            pans_i[i] = pans_i[h[iid]-1];
            goto ibl;
          }
          if(++iid >= nu) iid = 0;
        }
        h[iid] = i + 1; // need + 1 because for zero the while loop gives false..
        pans_i[i] = ++g;
        ibl:;
      }
    } else if(x_max == NA_INTEGER) { // fastver version based on range
      x_min -= 1;
      if(anyNA) {
        for (int i = 0, j; i != n; ++i) {
          j = (px[i] == NA_INTEGER) ? 0 : px[i]-x_min;
          if(h[j]) pans_i[i] = h[j];
          else pans_i[i] = h[j] = ++g;
        }
      } else {
        for (int i = 0, j; i != n; ++i) {
          j = px[i]-x_min;
          if(h[j]) pans_i[i] = h[j];
          else pans_i[i] = h[j] = ++g;
        }
      }
    } else {
      for (int i = 0; i != n; ++i) {
        id = HASH(px[i], K);
        while(h[id]) {
          if(px[h[id]-1] == px[i]) {
            pans_i[i] = pans_i[h[id]-1]; // h[id];
            goto ibbl;
          }
          if(++id >= M) id = 0;
        }
        h[id] = i + 1;
        pans_i[i] = ++g; // h[id];
        ibbl:;
      }
    }
  } break;
  case REALSXP: {
    const double *restrict px = REAL(x);
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      tpv.d = px[i] + 0.0; // R_IsNA(px[i]) ? NA_REAL : (R_IsNaN(px[i]) ? R_NaN : px[i]);
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto rbl;
        }
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      rbl:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *restrict px = COMPLEX(x);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp = px[i];
      if(C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r + 0.0;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i + 0.0;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto cbl;
        }
        if(++id >= M) id = 0; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      cbl:;
    }
  } break;
  case STRSXP: {
    const SEXP *restrict px = SEXPPTR(x);
    for (int i = 0; i != n; ++i) {
      id = HASH(((uintptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto sbl;
        }
        if(++id >= M) id = 0; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      sbl:;
    }
  } break;
  }
  R_Free(h);
  SEXP ngroups_sym = install("N.groups");
  setAttrib(ans_i, ngroups_sym, ScalarInteger(g));
  UNPROTECT(1);
  return ans_i;
}


SEXP dupVecIndexKeepNA(SEXP x) {
  const int n = length(x);
  int K = 0, tx = TYPEOF(x);
  size_t M;
  // if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
    bigint:;
    const size_t n2 = 2U * (size_t) n;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  } else if(tx == INTSXP) {
    if(isFactor(x) || inherits(x, "qG")) {
      tx = 1000;
      M = isFactor(x) ? (size_t)nlevels(x) + 2 : (size_t)asInteger(getAttrib(x, install("N.groups"))) + 2;
    } else {
      int *p = INTEGER(x);
      if(n > 10 && (NOGE(p[0], n) || NOGE(p[n/2], n) || NOGE(p[n-1], n))) goto bigint;
      M = (size_t)n;
    }
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  SEXP ans_i = PROTECT(allocVector(INTSXP, n));
  int *restrict pans_i = INTEGER(ans_i), g = 0;
  size_t id = 0;
  switch (tx) {
  case LGLSXP:
  case 1000: // This is for factors or logical vectors where the size of the table is known
  {
    const int *restrict px = INTEGER(x);
    for (int i = 0, j; i != n; ++i) {
      if(px[i] == NA_INTEGER) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      j = px[i];
      if(h[j]) pans_i[i] = h[j];
      else pans_i[i] = h[j] = ++g;
    }
  } break;
  case INTSXP: {
    const int *restrict px = INTEGER(x);
    if(M == (size_t)n) { // Faster version based on division hash...
      unsigned int iid = 0, nu = (unsigned)n;
      for (int i = 0; i != n; ++i) {
        if(px[i] == NA_INTEGER) {
          pans_i[i] = NA_INTEGER;
          continue;
        }
        iid = (unsigned)px[i];
        if(iid >= nu) iid %= nu; // iid = (px[i] < n) ? px[i] : px[i] % n; // HASH(px[i], K); // get the hash value of x[i]
        while(h[iid]) { // Check if this hash value has been seen before
          if(px[h[iid]-1] == px[i]) { // Get the element of x that produced his value. if x[i] is the same, assign it the same index.
            pans_i[i] = pans_i[h[iid]-1]; // h[id];
            goto ibl;
          } // else, we move forward to the next slot, until we find an empty one... We need to keep checking against the values,
          // because if we found the same value before, we would also have put it in another slot after the initial one with the same hash value.
          if(++iid >= nu) iid = 0; // ++iid; iid %= nu; // # nocov
        } // We put the index into the empty slot.
        h[iid] = i + 1; // need + 1 because for zero the while loop gives false..
        pans_i[i] = ++g; // h[id];
        ibl:;
      }
    } else {
      for (int i = 0; i != n; ++i) {
        if(px[i] == NA_INTEGER) {
          pans_i[i] = NA_INTEGER;
          continue;
        }
        id = HASH(px[i], K);
        while(h[id]) {
          if(px[h[id]-1] == px[i]) {
            pans_i[i] = pans_i[h[id]-1]; // h[id];
            goto ibbl;
          }
          if(++id >= M) id = 0; // ++id; id %= M;
        }
        h[id] = i + 1;
        pans_i[i] = ++g; // h[id];
        ibbl:;
      }
    }
  } break;
  case REALSXP: {
    const double *restrict px = REAL(x);
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      if(ISNAN(px[i])) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      tpv.d = px[i] + 0.0;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto rbl;
        }
        if(++id >= M) id = 0; // ++id; id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      rbl:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *restrict px = COMPLEX(x);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp = px[i];
      if(C_IsNA(tmp) || C_IsNaN(tmp)) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      tpv.d = tmp.r + 0.0;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i + 0.0;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto cbl;
        }
        if(++id >= M) id = 0; //++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      cbl:;
    }
  } break;
  case STRSXP: {
    const SEXP *restrict px = SEXPPTR(x);
    for (int i = 0; i != n; ++i) {
      if(px[i] == NA_STRING) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      id = HASH(((uintptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto sbl;
        }
        if(++id >= M) id = 0; //++id; id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g;
      sbl:;
    }
  } break;
  }
  R_Free(h);
  SEXP ngroups_sym = install("N.groups");
  setAttrib(ans_i, ngroups_sym, ScalarInteger(g));
  UNPROTECT(1);
  return ans_i;
}


// ****************************************
// Group Two Vectors in One Pass
// ****************************************
SEXP dupVecIndexTwoVectors(SEXP x, SEXP y) {

  int n = length(x), tx = TYPEOF(x), ty = TYPEOF(y), K, K2, anyNA = 0;
  if(length(y) != n) error("length of first two columns in the data must be the same");
  if(tx == CPLXSXP || ty == CPLXSXP) return R_NilValue;

  size_t M;
  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *restrict pans = INTEGER(ans);
  // Check if both are discrete
  int both_discr = (tx == LGLSXP || (tx == INTSXP && (isFactor(x) || inherits(x, "qG")))) &&
                   (ty == LGLSXP || (ty == INTSXP && (isFactor(y) || inherits(y, "qG"))));
  if(both_discr) {
    K = tx == LGLSXP ? 1 : isFactor(x) ? nlevels(x) : asInteger(getAttrib(x, install("N.groups")));
    K2 = ty == LGLSXP ? 1 : isFactor(y) ? nlevels(y) : asInteger(getAttrib(y, install("N.groups")));
    if(tx == LGLSXP || !inherits(x, "na.included")) {
      K += 1; anyNA += 1;
    }
    if(ty == LGLSXP || !inherits(y, "na.included")) {
      K2 += 1; anyNA += 1;
    }
    if((size_t)K * K2 <= (size_t)n * 3) {
      M = anyNA ? (size_t)(K+1) * (K2+1) + 1 : (size_t)K * K2 + 1; // + 1 because of zero indexing
    } else both_discr = 0;
  }

  if(!both_discr) {
    const size_t n2 = 2U * (size_t)n;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  }

  int *restrict h = (int*)R_Calloc(M, int), g = 0, hid = 0; // Table to save the hash values, table has size M
  size_t id = 0;

  if(both_discr) {
    const int *restrict px = INTEGER_RO(x), *restrict py = INTEGER_RO(y);
    if(anyNA == 0) {
      for (int i = 0; i != n; ++i) {
        id = px[i] + (py[i]-1) * K; // This assumes logical vectors (with 0 = FALSE) cannot have the "na.included" attribute
        if(h[id]) pans[i] = h[id];
        else pans[i] = h[id] = ++g;
      }
    } else {
      K += 1;
      for (int i = 0, xi, yi; i != n; ++i) {
        xi = (px[i] == NA_INTEGER) ? K : px[i]+1;
        yi = (py[i] == NA_INTEGER) ? K2 : py[i];
        id = xi + yi * K; // Problem: if logical xi = 0, yi = 1 and xi = K, yi = 0 give the same..
        if(h[id]) pans[i] = h[id];
        else pans[i] = h[id] = ++g;
      }
    }
  } else {
    if(tx == LGLSXP) tx = INTSXP;
    if(ty == LGLSXP) ty = INTSXP;

    // 6 cases: 3 same type and 3 different types
    if(tx == ty) { // same type
      switch(tx) {
        case INTSXP: {
          const int *restrict px = INTEGER_RO(x), *restrict py = INTEGER_RO(y);
          for (int i = 0; i != n; ++i) {
            id = HASH(px[i] + (64988430769U * py[i]), K); // Multiplication doesn't work here: too few unique values
            // Another large prime taken from https://oeis.org/wiki/Higher-order_prime_numbers
            while(h[id]) {
              hid = h[id]-1;
              if(px[hid] == px[i] && py[hid] == py[i]) {
                pans[i] = pans[hid];
                goto iibl;
              }
              if(++id >= M) id = 0;
            }
            h[id] = i + 1;
            pans[i] = ++g;
            iibl:;
          }
        } break;
        case STRSXP: {
          const SEXP *restrict px = SEXPPTR_RO(x), *restrict py = SEXPPTR_RO(y);
          for (int i = 0; i != n; ++i) {
            id = HASH(64988430769U * ((uintptr_t)px[i] & 0xffffffff) + ((uintptr_t)py[i] & 0xffffffff), K); // Best combination it seems
            while(h[id]) {
              hid = h[id]-1;
              if(px[hid] == px[i] && py[hid] == py[i]) {
                pans[i] = pans[hid];
                goto ssbl;
              }
              if(++id >= M) id = 0;
            }
            h[id] = i + 1;
            pans[i] = ++g;
            ssbl:;
          }
        } break;
        case REALSXP: {
          const double *restrict px = REAL_RO(x), *restrict py = REAL_RO(y);
          union uno tpx, tpy;
          // fill hash table with indices of 'table'
          for (int i = 0; i != n; ++i) {
            tpx.d = px[i] + 0.0; tpy.d = py[i] + 0.0;
            id = HASH((64988430769U * (tpx.u[0] + tpx.u[1])) + tpy.u[0] + tpy.u[1], K); // Best combination it seems
            while(h[id]) {
              hid = h[id]-1;
              if(REQUAL(px[hid], px[i]) && REQUAL(py[hid], py[i])) {
                pans[i] = pans[hid];
                goto rrbl;
              }
              if(++id >= M) id = 0;
            }
            h[id] = i + 1;
            pans[i] = ++g;
            rrbl:;
          }
        } break;
        default: error("Type %s is not supported.", type2char(tx));
      }
    } else { // different types
      // First case: integer and real
      if((tx == INTSXP && ty == REALSXP) || (tx == REALSXP && ty == INTSXP)) {
        const int *restrict pi = INTEGER_RO(tx == INTSXP ? x : y);
        const double *restrict pr = REAL_RO(tx == REALSXP ? x : y);
        union uno tpv;
        for (int i = 0; i != n; ++i) {
          tpv.d = pr[i] + 0.0;
          id = HASH((64988430769U * pi[i]) + tpv.u[0] + tpv.u[1], K); // Best combination it seems
          while(h[id]) {
            hid = h[id]-1;
            if(pi[hid] == pi[i] && REQUAL(pr[hid], pr[i])) {
              pans[i] = pans[hid];
              goto irbl;
            }
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          pans[i] = ++g;
          irbl:;
        }
        // Second case: real and string
      } else if ((tx == REALSXP && ty == STRSXP) || (tx == STRSXP && ty == REALSXP)) {
        const SEXP *restrict ps = SEXPPTR_RO(tx == STRSXP ? x : y);
        const double *restrict pr = REAL_RO(tx == REALSXP ? x : y);
        union uno tpv;
        for (int i = 0; i != n; ++i) {
          tpv.d = pr[i] + 0.0;
          id = HASH((tpv.u[0] + tpv.u[1]) * ((uintptr_t)ps[i] & 0xffffffff), K); // Best combination it seems
          while(h[id]) {
            hid = h[id]-1;
            if(REQUAL(pr[hid], pr[i]) && ps[hid] == ps[i]) { // Seems comparing reals is faster..
              pans[i] = pans[hid];
              goto srbl;
            }
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          pans[i] = ++g;
          srbl:;
        }
        // Third case: integer and string
      } else if((tx == INTSXP && ty == STRSXP) || (tx == STRSXP && ty == INTSXP)) {
        const int *restrict pi = INTEGER_RO(tx == INTSXP ? x : y);
        const SEXP *restrict ps = SEXPPTR_RO(tx == STRSXP ? x : y);
        for (int i = 0; i != n; ++i) {
          id = HASH(pi[i] * ((uintptr_t)ps[i] & 0xffffffff), K);
          while(h[id]) {
            hid = h[id]-1;
            if(pi[hid] == pi[i] && ps[hid] == ps[i]) {
              pans[i] = pans[hid];
              goto isbl;
            }
            if(++id >= M) id = 0;
          }
          h[id] = i + 1;
          pans[i] = ++g;
          isbl:;
        }
      } else error("Unsupported types: %s and %s", type2char(tx), type2char(ty));
    }
  }

  R_Free(h);
  SEXP ngroups_sym = install("N.groups");
  setAttrib(ans, ngroups_sym, ScalarInteger(g));
  UNPROTECT(1);
  return ans;
}


// TODO: Only one M calculation ?
// Think: If in the second grouping variable all entries are the same, you loop through the whole table for each value..

// TODO: Faster possibility indexing by grouping vector?? -> would need multiple hash tables through which complicates things,
// but could still end up being faster...

// Idea: instead of hasing index again, just distribute it fairly through multiplying with (M/ng)

// **************************************************
// This function adds a second vector to the grouping
// **************************************************
int dupVecSecond(int *restrict pidx, int *restrict pans_i, SEXP x, const int n, const int ng) {

  if(length(x) != n) error("Unequal length columns");
  int K = 0, tx = TYPEOF(x), anyNA = 1;
  size_t M;
  if (tx == INTSXP || tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
    if(tx == INTSXP && (isFactor(x) || inherits(x, "qG"))) {
      K = isFactor(x) ? nlevels(x)+1 : asInteger(getAttrib(x, install("N.groups")))+1;
      anyNA = !inherits(x, "na.included");
      if((size_t)K * ng <= (size_t)n * 3) {
        tx = 1000;
        M = (size_t)K * ng + 1;
      } else K = 0;
    }
    if(K == 0) {
      const size_t n2 = 2U * (size_t)n; // + ng
      M = 256;
      K = 8;
      while (M < n2) {
        M *= 2;
        K++;
      }
      // M += ng; // Here we add the number of previous groups...
    }
  } else if (tx == LGLSXP) {
    M = (size_t)ng * 3 + 1;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *restrict h = (int*)R_Calloc(M, int), g = 0, hid = 0; // Table to save the hash values, table has size M
  size_t id = 0;
  switch (tx) {
  case LGLSXP:
  {
    const int *restrict px = LOGICAL(x);
    for (int i = 0; i != n; ++i) {
      id = (px[i] == NA_LOGICAL) ? pidx[i] : pidx[i] + (px[i] + 1) * ng;
      if(h[id]) pans_i[i] = h[id];
      else pans_i[i] = h[id] = ++g;
    }
  } break;
  case 1000: // This is for factors if feasible...
  {
    const int *restrict px = INTEGER(x);
    if(anyNA) {
      for (int i = 0; i != n; ++i) {
        id = (px[i] == NA_INTEGER) ? pidx[i] : pidx[i] + px[i] * ng;
        if(h[id]) pans_i[i] = h[id];
        else pans_i[i] = h[id] = ++g;
      }
    } else {
      for (int i = 0; i != n; ++i) {
        id = pidx[i] + px[i] * ng;
        if(h[id]) pans_i[i] = h[id];
        else pans_i[i] = h[id] = ++g;
      }
    }
  } break;
  // TODO: Think further about this! Perhaps you can also do this totally differently with a second vector capturing the unique values of idx!
  // See again what Morgan does to his matrix of single groupings...

  // Note: In general, combining bitwise i.e. px[i] ^ pidx[i] seems slightly faster than multiplying (px[i] * pidx[i])...
  case INTSXP: {
    const int *restrict px = INTEGER(x);
    const unsigned int mult = (M-1) / ng; // -1 because C is zero indexed
    for (int i = 0; i != n; ++i) { // Check this... DATA group_by lon, lat main_cat, main_tag, main_tag_value: main_tag is the issue..
      id = (pidx[i]*mult) ^ HASH(px[i], K); // HASH((unsigned)px[i] * (unsigned)pidx[i], K) + pidx[i]; // Need multiplication here instead of bitwise, see your benchmark with 100 mio. obs where second group is just sample.int(1e4, 1e8, T), there bitwise is very slow!!
      while(h[id]) {  // However multiplication causes signed integer overflow... UBSAN error.
        hid = h[id]-1;
        if(pidx[hid] == pidx[i] && px[hid] == px[i]) { // Usually pidx has more distinct values...
          pans_i[i] = pans_i[hid];
          goto ibl;
        }
        if(++id >= M) id = 0; // ++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      ibl:;
    }
  } break;
  case REALSXP: {
    const double *restrict px = REAL(x);
    const unsigned int mult = (M-1) / ng; // -1 because C is zero indexed
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      tpv.d = px[i] + 0.0; // R_IsNA(px[i]) ? NA_REAL : (R_IsNaN(px[i]) ? R_NaN :px[i]);
      id = (pidx[i]*mult) ^ HASH(tpv.u[0] + tpv.u[1], K); // HASH((tpv.u[0] + tpv.u[1]) ^ pidx[i], K) + pidx[i]; // Note: This is much faster than just adding pidx[i] to the hash value...
      while(h[id]) { // Problem: This value might be seen before, but not in combination with that pidx value...
        hid = h[id]-1; // The issue here is that REQUAL(px[hid], px[i]) could be true but pidx[hid] == pidx[i] fails, although the same combination of px and pidx could be seen earlier before...
        if(pidx[hid] == pidx[i] && REQUAL(px[hid], px[i])) {
          pans_i[i] = pans_i[hid];
          goto rbl;
        }
        if(++id >= M) id = 0; //++id; id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      rbl:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *restrict px = COMPLEX(x);
    const unsigned int mult = (M-1) / ng; // -1 because C is zero indexed
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp = px[i];
      if(C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r + 0.0;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i + 0.0;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = (pidx[i]*mult) ^ HASH(u, K);
      while(h[id]) {
        hid = h[id]-1;
        if(pidx[hid] == pidx[i] && CEQUAL(px[hid], px[i])) {
          pans_i[i] = pans_i[hid];
          goto cbl;
        }
        if(++id >= M) id = 0; //++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      cbl:;
    }
  } break;
  case STRSXP: {
    const SEXP *restrict px = SEXPPTR(x);
    const unsigned int mult = (M-1) / ng; // -1 because C is zero indexed
    for (int i = 0; i != n; ++i) {
      id = (pidx[i]*mult) ^ HASH(((uintptr_t) px[i] & 0xffffffff), K); // HASH(((uintptr_t) px[i] & 0xffffffff) ^ pidx[i], K) + pidx[i];
      while(h[id]) {
        hid = h[id]-1;
        if(pidx[hid] == pidx[i] && px[hid] == px[i]) {
          pans_i[i] = pans_i[hid];
          goto sbl;
        }
        if(++id >= M) id = 0; //++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      sbl:;
    }
  } break;
  }
  R_Free(h);
  return g;
}


// ************************************************************************
// This function brings everything together for vectors or lists of vectors
// ************************************************************************
SEXP groupVec(SEXP X, SEXP starts, SEXP sizes) {

  int l = length(X), islist = TYPEOF(X) == VECSXP,
    start = asLogical(starts), size = asLogical(sizes), nprotect = 0;
  // Better not exceptions to fundamental algorithms, when a couple of user-level functions return qG objects...
  // if(islist == 0 && OBJECT(X) != 0 && inherits(X, "qG") && inherits(X, "na.included")) return X; // return "qG" objects
  const SEXP *px = islist ? SEXPPTR_RO(X) : &X;
  SEXP idx = islist == 0 ? dupVecIndex(X) : l > 1 ? dupVecIndexTwoVectors(px[0], px[1]) : dupVecIndex(px[0]);
  if(isNull(idx)) { // One of the vectors is complex valued
    idx = dupVecIndex(px[0]);
    l += 1; px -= 1;
  } else if(!(islist && l > 2) && start == 0 && size == 0) return idx; // l == 1 &&
  PROTECT(idx); ++nprotect;
  SEXP sym_ng = install("N.groups"), res;
  int ng = asInteger(getAttrib(idx, sym_ng)), n = length(idx);
  if(islist && l > 2) {
    SEXP ans = PROTECT(allocVector(INTSXP, n)); ++nprotect;
    int i = 2, *pidx = INTEGER(idx), *pans = INTEGER(ans);
    for( ; i < l; ++i) {
      if(ng == n) break;
      if(i % 2) {
        ng = dupVecSecond(pans, pidx, px[i], n, ng);
      } else {
        ng = dupVecSecond(pidx, pans, px[i], n, ng);
      }
    }
    res = i % 2 ? ans : idx;
    setAttrib(res, sym_ng, ScalarInteger(ng));
  } else res = idx;
  // Cumpoting group starts and sizes attributes
  if(start || size) {
    PROTECT(res); ++nprotect;
    int *pres = INTEGER(res);
    if(start && size) { // Protect res ??
      SEXP gs, st, starts_sym = install("starts"), sizes_sym = install("group.sizes");
      setAttrib(res, starts_sym, st = allocVector(INTSXP, ng));
      setAttrib(res, sizes_sym, gs = allocVector(INTSXP, ng));
      if(ng > 0) {
        int *pgs = INTEGER(gs), *pst = INTEGER(st);
        memset(pgs, 0, sizeof(int) * ng); --pgs;
        memset(pst, 0, sizeof(int) * ng); --pst;
        for(int i = 0; i != n; ++i) {
          ++pgs[pres[i]];
          if(pst[pres[i]] == 0) pst[pres[i]] = i + 1;
        }
      }
    } else if(start) {
      SEXP st, starts_sym = install("starts");
      setAttrib(res, starts_sym, st = allocVector(INTSXP, ng));
      if(ng > 0) {
        int *pst = INTEGER(st), k = 0;
        memset(pst, 0, sizeof(int) * ng); --pst;
        for(int i = 0; i != n; ++i) {
          if(pst[pres[i]] == 0) {
            pst[pres[i]] = i + 1;
            if(++k == ng) break;
          }
        }
      }
    } else {
      SEXP gs, sizes_sym = install("group.sizes");
      setAttrib(res, sizes_sym, gs = allocVector(INTSXP, ng));
      if(ng > 0) {
        int *pgs = INTEGER(gs);
        memset(pgs, 0, sizeof(int) * ng); --pgs;
        for(int i = 0; i != n; ++i) ++pgs[pres[i]];
      }
    }
  }
  UNPROTECT(nprotect);
  return res;
}


// This version is only for atomic vectors (factor generation)
SEXP groupAtVec(SEXP X, SEXP starts, SEXP naincl) {

  int start = asLogical(starts), nain = asLogical(naincl);
  // Note: These functions will give errors for unsupported types...
  SEXP idx = nain ? dupVecIndex(X) : dupVecIndexKeepNA(X);
  if(start == 0) return idx;
  PROTECT(idx);
  SEXP st, sym_ng = install("N.groups"), starts_sym = install("starts");
  int ng = asInteger(getAttrib(idx, sym_ng)), n = length(idx), *pidx = INTEGER(idx);
  setAttrib(idx, starts_sym, st = allocVector(INTSXP, ng));
  if(ng > 0) {
    int *pst = INTEGER(st), k = 0;
    memset(pst, 0, sizeof(int) * ng); --pst;
    if(nain) {
      for(int i = 0; i != n; ++i) {
        if(pst[pidx[i]] == 0) {
          pst[pidx[i]] = i + 1;
          if(++k == ng) break;
        }
      }
    } else {
      for(int i = 0; i != n; ++i) {
        if(pidx[i] != NA_INTEGER && pst[pidx[i]] == 0) {
          pst[pidx[i]] = i + 1;
          if(++k == ng) break;
        }
      }
    }
  }
  UNPROTECT(1);
  return idx;
}


// Same as dupVecIndex, but saves group starts and returns unique values
SEXP funiqueC(SEXP x) {
  const int n = length(x);
  if(n <= 1) return x;
  int K = 0, tx = TYPEOF(x);
  size_t M;
  // if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP) {
    bigint:;
    const size_t n2 = 2U * (size_t) n;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  } else if(tx == INTSXP) {
    if(isFactor(x) || inherits(x, "qG")) {
      tx = 1000;
      M = isFactor(x) ? (size_t)nlevels(x) + 2 : (size_t)asInteger(getAttrib(x, install("N.groups"))) + 2;
    } else {
      int *p = INTEGER(x);
      if(n > 10 && (NOGE(p[0], n) || NOGE(p[n/2], n) || NOGE(p[n-1], n))) goto bigint;
      M = (size_t)n;
    }
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  int *restrict st = (int*)R_alloc((tx == LGLSXP || tx == 1000) ? (int)M : n, sizeof(int));
  int g = 0, nprotect = 0;
  size_t id = 0;
  SEXP res = R_NilValue;
  switch (tx) {
  case LGLSXP:
  case 1000: // This is for factors or logical vectors where the size of the table is known
  {
    const int *restrict px = INTEGER(x);
    if(tx == 1000 && inherits(x, "na.included")) {
      for(int i = 0, k = (int)M-1, ng = k-1; i != n; ++i) {
        if(h[px[i]]) continue;
        h[px[i]] = 1;
        st[g] = i;
        if(++g == ng) break;
      }
    } else {
      int  ng = tx == LGLSXP ? 3 : (int)M-1;
      for(int i = 0, j, k = (int)M-1; i != n; ++i) {
        j = (px[i] == NA_INTEGER) ? k : px[i];
        if(h[j]) continue;
        h[j] = 1;
        st[g] = i;
        if(++g == ng) break;
      }
    }
    R_Free(h);
    if(g == n) return x;
    PROTECT(res = allocVector(tx == LGLSXP ? LGLSXP : INTSXP, g)); ++nprotect;
    int *restrict pres = INTEGER(res);
    for(int i = 0; i != g; ++i) pres[i] = px[st[i]];
  } break;
  case INTSXP: {
    const int *restrict px = INTEGER(x);
    if(M == (size_t)n) { // Faster version based on division hash...
      unsigned int iid = 0, nu = (unsigned)n;
      for (int i = 0; i != n; ++i) {
        iid = (unsigned)px[i];
        if(iid >= nu) iid %= nu;
        while(h[iid]) {
          if(px[h[iid]-1] == px[i]) goto ibl;
          if(++iid >= nu) iid = 0;
        }
        h[iid] = i + 1;
        st[g++] = i;
        ibl:;
      }
    } else {
      for (int i = 0; i != n; ++i) {
        id = HASH(px[i], K);
        while(h[id]) {
          if(px[h[id]-1] == px[i]) goto ibbl;
          if(++id >= M) id = 0;
        }
        h[id] = i + 1;
        st[g++] = i;
        ibbl:;
      }
    }
    R_Free(h);
    if(g == n) {
      UNPROTECT(nprotect);
      return x;
    }
    PROTECT(res = allocVector(INTSXP, g)); ++nprotect;
    int *restrict pres = INTEGER(res);
    for(int i = 0; i != g; ++i) pres[i] = px[st[i]];
  } break;
  case REALSXP: {
    const double *restrict px = REAL(x);
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      tpv.d = px[i] + 0.0;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[h[id]-1], px[i])) goto rbl;
        if(++id >= M) id = 0;
      }
      h[id] = i + 1;
      st[g++] = i;
      rbl:;
    }
    R_Free(h);
    if(g == n) {
      UNPROTECT(nprotect);
      return x;
    }
    PROTECT(res = allocVector(REALSXP, g)); ++nprotect;
    double *restrict pres = REAL(res);
    for(int i = 0; i != g; ++i) pres[i] = px[st[i]];
  } break;
  case CPLXSXP: {
    const Rcomplex *restrict px = COMPLEX(x);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp = px[i];
      if(C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r + 0.0;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i + 0.0;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(px[h[id]-1], px[i])) goto cbl;
        if(++id >= M) id = 0; // # nocov
      }
      h[id] = i + 1;
      st[g++] = i;
      cbl:;
    }
    R_Free(h);
    if(g == n) {
      UNPROTECT(nprotect);
      return x;
    }
    PROTECT(res = allocVector(CPLXSXP, g)); ++nprotect;
    Rcomplex *restrict pres = COMPLEX(res);
    for(int i = 0; i != g; ++i) pres[i] = px[st[i]];
  } break;
  case STRSXP: {
    const SEXP *restrict px = SEXPPTR(x);
    for (int i = 0; i != n; ++i) {
      id = HASH(((uintptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) goto sbl;
        if(++id >= M) id = 0; // # nocov
      }
      h[id] = i + 1;
      st[g++] = i;
      sbl:;
    }
    R_Free(h);
    if(g == n) {
      UNPROTECT(nprotect);
      return x;
    }
    PROTECT(res = allocVector(STRSXP, g)); ++nprotect;
    SEXP *restrict pres = SEXPPTR(res);
    for(int i = 0; i != g; ++i) pres[i] = px[st[i]];
  } break;
  }
  copyMostAttrib(x, res);
  UNPROTECT(nprotect);
  return res;
}

// TODO: fduplicated and any_duplicated: smart default methods...

// From the kit package...

/*
 *  Data.Frame
 */

// SEXP dupDataFrameR(SEXP x) { // move to matrix if possible
//
//   const SEXP *restrict px = SEXPPTR_RO(x);
//   const R_xlen_t len_x = xlength(x);
//   const R_xlen_t len_i = xlength(px[0]);
//   SEXP ans = R_NilValue;
//   SEXP mlv = PROTECT(allocMatrix(INTSXP, (int)len_i, (int)len_x));
//   for (R_xlen_t i = 0; i < len_x; ++i) {
//     memcpy(INTEGER(mlv)+i*len_i, INTEGER(PROTECT(dupVecIndexOnlyR(px[i]))), (unsigned)len_i*sizeof(int));
//   }
//   UNPROTECT((int)len_x);
//   const size_t n2 = 2U * (size_t) len_i;
//   size_t M = 256;
//   int K = 8;
//   while (M < n2) {
//     M *= 2;
//     K++;
//   }
//   R_xlen_t count = 0;
//   int *restrict h = (int*) R_Calloc(M, int);
//   const int *restrict v = INTEGER(mlv);
//   int *restrict pans = (int*) R_Calloc(len_i, int);
//   size_t id = 0;
//
//       for (R_xlen_t i = 0; i < len_i; ++i) {
//         R_xlen_t key = 0;
//         for (R_xlen_t j = 0; j < len_x; ++j) {
//           key ^= HASH(v[i+j*len_i],K)*97*(j+1);
//         }
//         id = HASH(key, K);
//         while (h[id]) {
//           for (R_xlen_t j = 0; j < len_x; ++j) {
//             if (v[h[id]-1+j*len_i] != v[i+j*len_i]) {
//               goto label1;
//             }
//           }
//           goto label2;
//           label1:;
//           id++; id %= M;
//         }
//         h[id] = (int) i + 1;
//         pans[i]++;
//         count++;
//         label2:;
//       }
//     R_Free(h);
//     UNPROTECT(1);
//     SEXP indx = PROTECT(allocVector(INTSXP, count));
//     int ct = 0;
//     int *restrict py = INTEGER(indx);
//     for (int i = 0; ct < count; ++i) {
//       if (pans[i]) {
//         py[ct++] = i;
//       }
//     }
//     SEXP output = PROTECT(subSetRowDataFrame(x, indx));
//     R_Free(pans);
//     UNPROTECT(2);
//     return output;
// }


/*
 *  Data.Frame
 */

// SEXP dupLenDataFrameR(SEXP x) {
//   const SEXP *restrict px = SEXPPTR_RO(x);
//   const R_xlen_t len_x = xlength(x);
//   // bool allT = true;
//   // const SEXPTYPE t0 = UTYPEOF(px[0]);
//   // for (int i = 1; i < len_x; ++i) {
//   //   if (UTYPEOF(px[i]) != t0) {
//   //     allT = false;
//   //     break;
//   //   }
//   // }
//   // if (allT) {
//   //   SEXP output = PROTECT(dupLenMatrixR(PROTECT(dfToMatrix(x))));
//   //   UNPROTECT(2);
//   //   return output;
//   // }
//   const R_xlen_t len_i = xlength(px[0]);
//   SEXP mlv = PROTECT(allocMatrix(INTSXP, (int)len_i, (int)len_x));
//   for (R_xlen_t i = 0; i < len_x; ++i) {
//     memcpy(INTEGER(mlv)+i*len_i, INTEGER(PROTECT(dupVecIndexOnlyR(px[i], ScalarLogical(false)))), (unsigned)len_i*sizeof(int));
//   }
//   UNPROTECT((int)len_x);
//   const size_t n2 = 2U * (size_t) len_i;
//   size_t M = 256;
//   int K = 8;
//   while (M < n2) {
//     M *= 2;
//     K++;
//   }
//   R_xlen_t count = 0;
//   int *restrict h = (int*) R_Calloc(M, int);
//   const int *restrict v = INTEGER(mlv);
//   size_t id = 0;
//   for (R_xlen_t i = 0; i < len_i; ++i) {
//     R_xlen_t key = 0;
//     for (R_xlen_t j = 0; j < len_x; ++j) {
//       key ^= HASH(v[i+j*len_i],K)*97*(j+1);
//     }
//     id = HASH(key, K);
//     while (h[id]) {
//       for (R_xlen_t j = 0; j < len_x; ++j) {
//         if (v[h[id]-1+j*len_i] != v[i+j*len_i]) {
//           goto label1;
//         }
//       }
//       goto label2;
//       label1:;
//       id++; id %= M;
//     }
//     h[id] = (int) i + 1;
//     count++;
//     label2:;
//   }
//   R_Free(h);
//   UNPROTECT(1);
//   return ScalarInteger(count);
// }
