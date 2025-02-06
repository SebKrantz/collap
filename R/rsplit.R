
# fsplit <- function(x, f, drop, ...) if(drop && is.factor(f))
#   split(x, .Call(Cpp_fdroplevels, f, !inherits(f, "na.included")), drop = FALSE, ...) else
#     split(x, qF(f), drop = FALSE, ...)

t_list2 <- function(x) .Call(Cpp_mctl, do.call(rbind, x), TRUE, 0L)

# This is for export
t_list <- function(l) {
  lmat <- do.call(rbind, l)
  dn <- dimnames(lmat)
  res <- .Call(Cpp_mctl, lmat, !is.null(dn[[2L]]), 0L)
  if(length(rn <- dn[[1L]])) res <- lapply(res, `names<-`, rn)
  .Call(C_copyMostAttrib, res, l)
}


rsplit <- function(x, ...) UseMethod("rsplit")

rsplit.default <- function(x, fl, drop = TRUE, flatten = FALSE, use.names = TRUE, ...) { # , check = TRUE
  # if(is.matrix(x) && !inherits(x, "matrix")) return(rsplit.matrix(x, fl, drop, flatten, use.names, ...))
  if(is.atomic(fl) || flatten || is_GRP(fl)) return(gsplit(x, fl, use.names, drop = drop, ...))
  attributes(fl) <- NULL
  # if(check) fl <- lapply(fl, qF) # necessary ? -> split.default is actually faster on non-factor variables !
  rspl <- function(y, fly) {
    if(length(fly) == 1L) return(gsplit(y, fly[[1L]], use.names, drop = drop, ...))
    mapply(rspl, y = gsplit(y, fly[[1L]], use.names, drop = drop, ...),
           fly = t_list2(lapply(fly[-1L], gsplit, fly[[1L]], use.names, drop = drop, ...)), SIMPLIFY = FALSE) # Possibility to avoid transpose ? C_subsetDT ??
  }
  rspl(x, fl)
}

# Matrix method: requested in https://github.com/ycroissant/plm/issues/33
split_mat <- function(x, fl, dd, ...) {
  ssfun <- if(dd) function(i) x[i, , drop = TRUE] else function(i) x[i, , drop = FALSE]
  lapply(gsplit(NULL, fl, ...), ssfun)
}

rsplit.matrix <- function(x, fl, drop = TRUE, flatten = FALSE, use.names = TRUE, drop.dim = FALSE, ...) {
  if(is.atomic(fl) || flatten || is_GRP(fl)) return(split_mat(x, fl, drop.dim, use.names, drop = drop, ...))
  attributes(fl) <- NULL
  rspl <- function(y, fly) {
    if(length(fly) == 1L) return(split_mat(y, fly[[1L]], drop.dim, use.names, drop = drop, ...))
    mapply(rspl, y = split_mat(y, fly[[1L]], drop.dim, use.names, drop = drop, ...),
           fly = t_list2(lapply(fly[-1L], gsplit, fly[[1L]], use.names, drop = drop, ...)), SIMPLIFY = FALSE)
  }
  rspl(x, fl)
}

rsplit.zoo <- function(x, ...) if(is.matrix(x)) rsplit.matrix(x, ...) else rsplit.default(x, ...)
rsplit.units <- rsplit.zoo

# From stackoverflow package:
# rsplit <- function (x, by, drop = FALSE)
# {
#   if (is.atomic(by))
#     return(split(x, by, drop = drop))
#   attributes(by) <- NULL
#   if (length(by) == 1L)
#     return(split(x, by[[1L]], drop = drop))
#   mapply(rsplit, x = split(x, by[[1L]], drop = drop), by = t(lapply(by[-1L], split, by[[1L]], drop = drop)), drop = drop,
#          SIMPLIFY = FALSE)
# }

rsplit.data.frame <- function(x, by, drop = TRUE, flatten = FALSE, # check = TRUE,
                              cols = NULL, keep.by = FALSE, simplify = TRUE,
                              use.names = TRUE, nest = FALSE,  ...) {

  if(is.call(by)) {
    nam <- attr(x, "names")
    if(length(by) == 3L) {
      byn <- ckmatch(all.vars(by[[3L]]), nam)
      cols <- ckmatch(all.vars(by[[2L]]), nam)
    } else { # keep.by always added: Same behavior as L or W !!
      byn <- ckmatch(all.vars(by), nam)
      if(!(is.null(cols) && keep.by))
        cols <- if(is.null(cols)) -byn else cols2int(cols, x, nam, FALSE)
    }
    by <- .subset(x, byn)
    if(length(cols)) x <- fcolsubset(x, if(keep.by) c(byn, cols) else cols, TRUE)
  } else if(length(cols))
    x <- fcolsubset(x, cols2int(cols, x, attr(x, "names"), FALSE), TRUE)

  if(simplify && length(unclass(x)) == 1L)
    return(rsplit.default(.subset2(x, 1L), by, drop, flatten, use.names, ...))  # , check
  # Note there is a data.table method: split.data.table, which can also do recursive splitting..

  j <- seq_along(unclass(x))
  rn <- attr(x, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") {
    gsplit_DF <- function(x, f, ...)
      lapply(gsplit(NULL, f, use.names, drop = drop, ...),
             function(i) .Call(C_subsetDT, x, i, j, FALSE)) # .Call, .NAME = C_subsetDT, j, FALSE) -> doesn't work!
  } else {
    gsplit_DF <- function(x, f, ...) {
      rown <- attr(x, "row.names") # Need to do this, handing down from the function body doesn't work
      lapply(gsplit(NULL, f, use.names, drop = drop, ...),
             function(i) `attr<-`(.Call(C_subsetDT, x, i, j, FALSE), "row.names", rown[i]))
    }
  }

  if(nest) return(ftransform(unique(data.frame(by)), data = gsplit_DF(x, by)))

  if(is.atomic(by) || flatten || is_GRP(by)) return(gsplit_DF(x, by, ...))

  attributes(by) <- NULL
  # if(check) by <- lapply(by, qF) # necessary ?
  rspl_DF <- function(y, fly) {
    if(length(fly) == 1L) return(gsplit_DF(y, fly[[1L]], ...))
    mapply(rspl_DF, y = gsplit_DF(y, fly[[1L]], ...),
           fly = t_list2(lapply(fly[-1L], gsplit, fly[[1L]], use.names, drop = drop, ...)), SIMPLIFY = FALSE) # Possibility to avoid transpose ?
  }                # use C_subsetDT here as well ??? what is faster ???
  rspl_DF(x, by)
}

