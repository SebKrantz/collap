#include "collapse_c.h"

void fcumsum_double_impl(double *pout, double *px, int ng, int *pg, int narm, int fill, int l) {
  if(ng == 0) {
    if(narm <= 0) {
      pout[0] = px[0];
      for(int i = 1; i != l; ++i) pout[i] = pout[i-1] + px[i];
    } else if(fill) {
      pout[0] = ISNAN(px[0]) ? 0.0 : px[0];
      for(int i = 1; i != l; ++i) pout[i] = pout[i-1] + (ISNAN(px[i]) ? 0.0 : px[i]);
    } else {
      double last = 0;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) pout[i] = px[i];
        else pout[i] = last += px[i];
      }
    }
  } else {
    double last[ng+1]; // Also pass pointer to function ??
    memset(last, 0.0, sizeof(double) * (ng+1));
    if(narm <= 0) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + px[i];
    } else if(fill) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + (ISNAN(px[i]) ? 0.0 : px[i]);
    } else {
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) pout[i] = px[i];
        else last[pg[i]] = pout[i] = last[pg[i]] + px[i];
      }
    }
  }
}

void fcumsum_double_impl_order(double *pout, double *px, int ng, int *pg, int *po, int narm, int fill, int l) {
  if(ng == 0) {
    if(narm <= 0) {
      --pout; --px;
      pout[po[0]] = px[po[0]];
      for(int i = 1; i != l; ++i) pout[po[i]] = pout[po[i-1]] + px[po[i]];
    } else if(fill) {
      --pout; --px;
      pout[po[0]] = ISNAN(px[po[0]]) ? 0.0 : px[po[0]];
      for(int i = 1; i != l; ++i) pout[po[i]] = pout[po[i-1]] + (ISNAN(px[po[i]]) ? 0.0 : px[po[i]]);
    } else {
      double last = 0;
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        if(ISNAN(px[poi])) pout[poi] = px[poi];
        else pout[poi] = last += px[poi];
      }
    }
  } else {
    double last[ng+1]; // Also pass pointer to function ??
    memset(last, 0.0, sizeof(double) * (ng+1));
    if(narm <= 0) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    } else if(fill) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        last[pg[poi]] = pout[poi] = last[pg[poi]] + (ISNAN(px[poi]) ? 0.0 : px[poi]);
      }
    } else {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        if(ISNAN(px[poi])) pout[poi] = px[poi];
        else last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    }
  }
}


void fcumsum_int_impl(int *pout, int *px, int ng, int *pg, int narm, int fill, int l) {
  if(ng == 0) {
    double ckof;
    if(narm <= 0) {
      int i = 1;
      ckof = pout[0] = px[0];
      if(px[0] == NA_INTEGER) --i;
      for( ; i != l; ++i) {
        if(px[i] == NA_INTEGER) {
          for( ; i != l; ++i) pout[i] = NA_INTEGER;
          break;
        }
        pout[i] = ckof += px[i];
      }
    } else if(fill) {
      ckof = pout[0] = (px[0] == NA_INTEGER) ? 0 : px[0];
      for(int i = 1; i != l; ++i) {
        if(px[i] != NA_INTEGER) ckof += px[i];
        pout[i] = ckof;
      }
    } else {
      ckof = 0.0;
      for(int i = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) pout[i] = NA_INTEGER;
        else pout[i] = ckof += px[i];
      }
    }
    if(ckof > INT_MAX || ckof < INT_MIN || (narm && ckof == NA_INTEGER))
      error("Integer overflow, use fcumsum(as.numeric(x))");
  } else {
    int last[ng+1]; // Also pass pointer to function ??
    memset(last, 0, sizeof(int) * (ng+1));
    if(narm <= 0) {
      for(int i = 0, lsi; i != l; ++i) {
        if(px[i] == NA_INTEGER) {
          pout[i] = last[pg[i]] = NA_INTEGER;
          continue;
        }
        lsi = last[pg[i]];
        if(lsi == NA_INTEGER) pout[i] = NA_INTEGER;
        else last[pg[i]] = pout[i] = lsi + px[i];
      }
    } else if(fill) {
      for(int i = 0; i != l; ++i) last[pg[i]] = pout[i] = last[pg[i]] + (px[i] == NA_INTEGER ? 0 : px[i]);
    } else {
      for(int i = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) pout[i] = NA_INTEGER;
        else last[pg[i]] = pout[i] = last[pg[i]] + px[i];
      }
    }
  }
}

void fcumsum_int_impl_order(int *pout, int *px, int ng, int *pg, int *po, int narm, int fill, int l) {
  if(ng == 0) {
    double ckof;
    if(narm <= 0) {
      int i = 1, poi;
      ckof = pout[po[0]-1] = px[po[0]-1];
      if(px[po[0]-1] == NA_INTEGER) --i;
      for( ; i != l; ++i) {
        poi = po[i]-1;
        if(px[poi] == NA_INTEGER) {
          for( ; i != l; ++i) pout[po[i]-1] = NA_INTEGER;
          break;
        }
        pout[poi] = ckof += px[poi];
      }
    } else if(fill) {
      ckof = pout[po[0]-1] = (px[po[0]-1] == NA_INTEGER) ? 0 : px[po[0]-1];
      for(int i = 1, poi; i != l; ++i) {
        poi = po[i]-1;
        if(px[poi] != NA_INTEGER) ckof += px[poi];
        pout[poi] = ckof;
      }
    } else {
      ckof = 0.0;
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        if(px[poi] == NA_INTEGER) pout[poi] = NA_INTEGER;
        else pout[poi] = ckof += px[poi];
      }
    }
    if(ckof > INT_MAX || ckof < INT_MIN || (narm && ckof == NA_INTEGER))
      error("Integer overflow, use fcumsum(as.numeric(x))");
  } else {
    int last[ng+1]; // Also pass pointer to function ??
    memset(last, 0, sizeof(int) * (ng+1));
    if(narm <= 0) {
      for(int i = 0, poi, lsi; i != l; ++i) {
        poi = po[i]-1;
        if(px[poi] == NA_INTEGER) {
          pout[poi] = last[pg[poi]] = NA_INTEGER;
          continue;
        }
        lsi = last[pg[poi]];
        if(lsi == NA_INTEGER) pout[poi] = NA_INTEGER;
        else last[pg[poi]] = pout[poi] = lsi + px[poi];
      }
    } else if(fill) {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        last[pg[poi]] = pout[poi] = last[pg[poi]] + (px[poi] == NA_INTEGER ? 0 : px[poi]);
      }
    } else {
      for(int i = 0, poi; i != l; ++i) {
        poi = po[i]-1;
        if(px[poi] == NA_INTEGER) pout[poi] = NA_INTEGER;
        else last[pg[poi]] = pout[poi] = last[pg[poi]] + px[poi];
      }
    }
  }
}

SEXP fcumsumC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), fill = asLogical(Rfill), *pg = INTEGER(g),
    ord  = length(o) > 1, *po = ord ? INTEGER(o) : pg;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng > 0 && l != length(g)) error("length(g) must match length(x)");
  if(ord && l != length(o)) error("length(o) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(tx, l));
  switch(tx) {
  case REALSXP:
    if(ord) fcumsum_double_impl_order(REAL(out), REAL(x), ng, pg, po, narm, fill, l);
    else fcumsum_double_impl(REAL(out), REAL(x), ng, pg, narm, fill, l);
    break;
  case INTSXP:
    if(ord) fcumsum_int_impl_order(INTEGER(out), INTEGER(x), ng, pg, po, narm, fill, l);
    else fcumsum_int_impl(INTEGER(out), INTEGER(x), ng, pg, narm, fill, l);
    break;
  default: error("Unsupported SEXP type");
  }
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}

SEXP fcumsummC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
     ng = asInteger(Rng), narm = asLogical(Rnarm), fill = asLogical(Rfill), *pg = INTEGER(g),
     ord  = length(o) > 1, *po = ord ? INTEGER(o) : pg;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng > 0 && l != length(g)) error("length(g) must match nrow(x)");
  if(ord && l != length(o)) error("length(o) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(tx, l * col));
  switch(tx) {
  case REALSXP: {
    double *px = REAL(x), *pout = REAL(out);
    if(ord) for(int j = 0; j != col; ++j) fcumsum_double_impl_order(pout + j*l, px + j*l, ng, pg, po, narm, fill, l);
    else for(int j = 0; j != col; ++j) fcumsum_double_impl(pout + j*l, px + j*l, ng, pg, narm, fill, l);
    break;
  }
  case INTSXP: {
    int *px = INTEGER(x), *pout = INTEGER(out);
    if(ord) for(int j = 0; j != col; ++j) fcumsum_int_impl_order(pout + j*l, px + j*l, ng, pg, po, narm, fill, l);
    else for(int j = 0; j != col; ++j) fcumsum_int_impl(pout + j*l, px + j*l, ng, pg, narm, fill, l);
    break;
  }
  default: error("Unsupported SEXP type");
  }
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}

SEXP fcumsumlC(SEXP x, SEXP Rng, SEXP g, SEXP o, SEXP Rnarm, SEXP Rfill) {
  int l = length(x);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
  for(int j = 0; j != l; ++j) pout[j] = fcumsumC(px[j], Rng, g, o, Rnarm, Rfill);
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}
