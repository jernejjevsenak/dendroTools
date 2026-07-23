#' boot_fast_prep
#'
#' Prepare the shared objects for the vectorised bootstrap. Called once per
#' search, before any window is evaluated.
#'
#' The bootstrap is rewritten as matrix algebra. If C[i, r] counts how many
#' times observation i was drawn into resample r, then every sum required by
#' a correlation is a weighted sum over the original observations, and the
#' whole set of replicates becomes one matrix product. Quantities that depend
#' only on the response are computed here and reused for every window.
#'
#' @param y numeric vector, the response (chronology) values
#' @param boot_n number of bootstrap replicates
#' @param cor_method "pearson", "spearman" or "kendall"
#' @param cor_na_use the missing-value rule, as passed to stats::cor()
#'
#' @return a list used by boot_fast_replicates()
#' @keywords internal
#'
boot_fast_prep <- function(y, boot_n, cor_method, cor_na_use = "everything") {

  n <- length(y)

  # resample indices are drawn once and shared by every window, so that
  # windows are compared on the same resampled years
  idx <- matrix(sample.int(n, n * boot_n, replace = TRUE), nrow = n)

  Cmat <- matrix(0, n, boot_n)
  for (r in seq_len(boot_n)) Cmat[, r] <- tabulate(idx[, r], nbins = n)

  p <- list(Cmat = Cmat, y = y, n = n, boot_n = boot_n,
            cor_method = cor_method, cor_na_use = cor_na_use,
            y_ties = anyDuplicated(y[!is.na(y)]) > 0)

  if (cor_method == "kendall") {

    s <- sign(outer(y, y, "-"))
    s[is.na(s)] <- 0
    p$SGNy <- s

  } else if (cor_method == "spearman") {

    a <- outer(y, y, ">")
    a[is.na(a)] <- FALSE
    p$ALTy <- a * 1.0
  }

  p
}


#' boot_fast_replicates
#'
#' All bootstrap replicates of the correlation between the response and one
#' aggregated climate window, computed without resampling loops.
#'
#' Missing values follow the rule recorded in the prepared object. For the
#' complete-case rules this is deletion within each replicate, which for two
#' variables is identical to pairwise deletion. For "everything" and "all.obs"
#' a replicate that draws any incomplete observation returns NA, matching
#' stats::cor().
#'
#' @param p object returned by boot_fast_prep()
#' @param x numeric vector, one aggregated climate window
#' @param min_complete minimum number of complete pairs required in a replicate
#'
#' @return numeric vector of length boot_n
#' @keywords internal
#'
boot_fast_replicates <- function(p, x, min_complete = 4L) {

  y <- p$y
  keep <- !is.na(x) & !is.na(y)
  full <- all(keep)

  Ce <- if (full) p$Cmat else { z <- p$Cmat; z[!keep, ] <- 0; z }
  nn <- if (full) rep(p$n, p$boot_n) else colSums(Ce)

  # ---------------- Pearson ----------------
  if (p$cor_method == "pearson") {

    # centring does not change a correlation but keeps the sums small
    # relative to their differences, which avoids cancellation error in
    # windows with little variance
    xc <- x - mean(x[keep])
    yc <- y - mean(y[keep])

    x0 <- xc; x0[!keep] <- 0
    y0 <- yc; y0[!keep] <- 0

    # one matrix product yields all six sums for all replicates
    B <- cbind(as.numeric(keep), y0, y0 * y0, x0, x0 * x0, x0 * y0)
    S <- crossprod(Ce, B)

    nn  <- S[, 1]
    sy  <- S[, 2]; syy <- S[, 3]
    sx  <- S[, 4]; sxx <- S[, 5]
    sxy <- S[, 6]

    out <- (nn * sxy - sx * sy) /
      sqrt((nn * sxx - sx^2) * (nn * syy - sy^2))

  # ---------------- Spearman ----------------
  } else if (p$cor_method == "spearman") {

    ALTy <- p$ALTy
    if (!full) ALTy[!keep, ] <- 0

    xx <- x; xx[!keep] <- NA
    ALTx <- outer(xx, xx, ">")
    ALTx[is.na(ALTx)] <- FALSE
    ALTx <- ALTx * 1.0

    x_ties <- anyDuplicated(x[keep]) > 0

    if (!p$y_ties && !x_ties) {

      # with untied originals, the only tied values in a resample are the
      # repeated draws of the same observation
      Ry <- (ALTy %*% Ce) + (Ce + 1) / 2
      Rx <- (ALTx %*% Ce) + (Ce + 1) / 2

    } else {

      eq_counts <- function(v) {
        vv <- v; vv[!keep] <- NA
        a <- outer(vv, vv, "=="); a[is.na(a)] <- FALSE
        (a * 1.0) %*% Ce
      }
      Ry <- (ALTy %*% Ce) + (eq_counts(y) + 1) / 2
      Rx <- (ALTx %*% Ce) + (eq_counts(x) + 1) / 2
    }

    sy  <- colSums(Ce * Ry);      sx  <- colSums(Ce * Rx)
    syy <- colSums(Ce * Ry * Ry); sxx <- colSums(Ce * Rx * Rx)
    sxy <- colSums(Ce * Ry * Rx)

    out <- (nn * sxy - sx * sy) /
      sqrt((nn * sxx - sx^2) * (nn * syy - sy^2))

  # ---------------- Kendall (tau-b) ----------------
  } else if (p$cor_method == "kendall") {

    SGNy <- p$SGNy
    if (!full) SGNy[!keep, ] <- 0

    xx <- x; xx[!keep] <- NA
    SGNx <- sign(outer(xx, xx, "-"))
    SGNx[is.na(SGNx)] <- 0

    # concordant minus discordant pairs is a quadratic form in the counts
    Q <- SGNy * SGNx
    S <- 0.5 * colSums(Ce * (Q %*% Ce))

    n0 <- nn * (nn - 1) / 2

    x_ties <- anyDuplicated(x[keep]) > 0

    if (!p$y_ties && !x_ties) {

      nt <- colSums(Ce * (Ce - 1)) / 2
      n1 <- nt; n2 <- nt

    } else {

      tie_term <- function(v) {
        vv <- v; vv[!keep] <- NA
        a <- outer(vv, vv, "=="); a[is.na(a)] <- FALSE
        G <- (a * 1.0) %*% Ce
        colSums(Ce * (G - 1) / 2)
      }
      n1 <- tie_term(y); n2 <- tie_term(x)
    }

    out <- S / sqrt((n0 - n1) * (n0 - n2))

  } else {

    stop("cor_method must be 'pearson', 'spearman' or 'kendall'.")
  }

  out[nn < min_complete] <- NA_real_
  out[!is.finite(out)] <- NA_real_

  # stats::cor() with use = "everything" (or "all.obs") returns NA as soon as
  # a single missing value is present, so any replicate that happened to draw
  # an incomplete observation must be NA. The complete-case rules
  # ("complete.obs", "na.or.complete", "pairwise.complete.obs") are already
  # satisfied by the masking above, which for two variables is identical to
  # pairwise deletion.
  if (!full && p$cor_na_use %in% c("everything", "all.obs")) {
    drawn_na <- as.numeric(crossprod(p$Cmat, as.numeric(!keep)))
    out[drawn_na > 0] <- NA_real_
  }

  out
}


#' boot_fast_norm_inter
#'
#' Interpolated order statistic used by boot::boot.ci() for the percentile and
#' basic intervals. Reimplemented because boot:::norm.inter is not exported.
#'
#' @param t numeric vector of bootstrap replicates
#' @param alpha numeric vector of probabilities
#'
#' @return numeric vector of interpolated quantiles
#' @keywords internal
#'
boot_fast_norm_inter <- function(t, alpha) {

  t <- t[is.finite(t)]
  R <- length(t)
  rk <- (R + 1) * alpha
  k <- trunc(rk)

  out <- numeric(length(k))
  tstar <- sort(t)

  ints <- (k == rk)
  if (any(ints)) out[ints] <- tstar[k[ints]]
  out[k == 0] <- tstar[1L]
  out[k == R] <- tstar[R]

  temp <- which(!ints & k != 0 & k != R)
  if (length(temp)) {
    a1 <- stats::qnorm(alpha[temp])
    a2 <- stats::qnorm(k[temp] / (R + 1))
    a3 <- stats::qnorm((k[temp] + 1) / (R + 1))
    tk  <- tstar[k[temp]]
    tk1 <- tstar[k[temp] + 1L]
    out[temp] <- tk + (a1 - a2) / (a3 - a2) * (tk1 - tk)
  }

  out
}


#' boot_fast_ci
#'
#' Confidence interval from bootstrap replicates, reproducing boot::boot.ci()
#' for the "norm", "perc" and "basic" interval types.
#'
#' @param t numeric vector of bootstrap replicates
#' @param t0 the statistic computed on the original sample
#' @param conf confidence level
#' @param type "norm", "perc" or "basic"
#'
#' @return numeric vector of length two, lower and upper limit
#' @keywords internal
#'
boot_fast_ci <- function(t, t0, conf = 0.95, type = "norm") {

  t <- t[is.finite(t)]

  if (length(t) < 2 || !is.finite(t0)) return(c(NA_real_, NA_real_))

  if (type == "norm") {

    bias <- mean(t) - t0
    merr <- stats::sd(t) * stats::qnorm((1 + conf) / 2)
    c(t0 - bias - merr, t0 - bias + merr)

  } else if (type == "perc") {

    qq <- boot_fast_norm_inter(t, (1 + c(-conf, conf)) / 2)
    c(qq[1], qq[2])

  } else if (type == "basic") {

    qq <- boot_fast_norm_inter(t, (1 + c(conf, -conf)) / 2)
    c(2 * t0 - qq[1], 2 * t0 - qq[2])

  } else {

    c(NA_real_, NA_real_)
  }
}


#' boot_fast_reprep
#'
#' Rebuild a prepared object around a different response vector while reusing
#' the resample count matrix. Needed for the primary-versus-control
#' correlation in a partial correlation, which does not involve the response.
#'
#' @param p object returned by boot_fast_prep()
#' @param y_new numeric vector to use as the response side
#'
#' @return a prepared object with the same resamples
#' @keywords internal
#'
boot_fast_reprep <- function(p, y_new) {

  q <- list(Cmat = p$Cmat, y = y_new, n = p$n, boot_n = p$boot_n,
            cor_method = p$cor_method, cor_na_use = p$cor_na_use,
            y_ties = anyDuplicated(y_new[!is.na(y_new)]) > 0)

  if (p$cor_method == "kendall") {

    s <- sign(outer(y_new, y_new, "-"))
    s[is.na(s)] <- 0
    q$SGNy <- s

  } else if (p$cor_method == "spearman") {

    a <- outer(y_new, y_new, ">")
    a[is.na(a)] <- FALSE
    q$ALTy <- a * 1.0
  }

  q
}


#' boot_fast_pcor_replicates
#'
#' All bootstrap replicates of the partial correlation between the response
#' and a primary climate window, controlling for the paired control window.
#'
#' Uses the identity
#'   r_xy.z = (r_xy - r_xz * r_yz) / (sqrt(1 - r_xz^2) * sqrt(1 - r_yz^2))
#' which reproduces psych::partial.r() exactly, so the three pairwise
#' correlations can each be obtained from the vectorised bootstrap.
#'
#' @param p object returned by boot_fast_prep(), prepared on the response
#' @param x_primary numeric vector, the aggregated primary climate window
#' @param x_control numeric vector, the paired aggregated control window
#' @param min_complete minimum number of complete cases required
#'
#' @return numeric vector of length boot_n
#' @keywords internal
#'
boot_fast_pcor_replicates <- function(p, x_primary, x_control,
                                      min_complete = 4L) {

  y <- p$y

  # Under the complete-case rules every pairwise correlation must be built
  # from the same rows, so the three series are masked jointly first. Under
  # "pairwise.complete.obs" each pair keeps its own complete cases, which is
  # what psych::partial.r() does when it forms the correlation matrix.
  if (p$cor_na_use %in% c("complete.obs", "na.or.complete")) {

    keep3 <- !is.na(y) & !is.na(x_primary) & !is.na(x_control)
    y[!keep3] <- NA
    x_primary[!keep3] <- NA
    x_control[!keep3] <- NA
    p <- boot_fast_reprep(p, y)
  }

  r_xy <- boot_fast_replicates(p, x_primary, min_complete = min_complete)
  r_xz <- boot_fast_replicates(p, x_control, min_complete = min_complete)

  # the primary-control correlation does not involve the response, so the
  # response-side precomputation cannot be reused here
  p_yz <- boot_fast_reprep(p, x_primary)
  r_yz <- boot_fast_replicates(p_yz, x_control, min_complete = min_complete)

  out <- (r_xy - r_xz * r_yz) /
    (sqrt(1 - r_xz^2) * sqrt(1 - r_yz^2))

  out[!is.finite(out)] <- NA_real_
  out
}
