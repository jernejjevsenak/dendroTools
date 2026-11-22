#' Bootstrap statistic for both partial correlations of a seascorr window
#'
#' Internal helper of the seascorr functions with pcor_control = TRUE. Returns
#' the partial correlation of the response (x) with the primary window (y),
#' controlling for the control window (z), and the partial correlation of the
#' response with the control window, controlling for the primary window. The
#' first value is computed exactly as boot_f_pcor() computes it, so one
#' bootstrap run serves both surfaces and the primary results do not change.
#'
#' @keywords internal
#' @noRd

boot_f_pcor_both <- function(data, indices, cor.type,
                             pcor.use = "pairwise.complete.obs") {

  d <- data[indices, , drop = FALSE]

  if (sum(stats::complete.cases(d[, c("x", "y", "z")])) < 4) {
    return(c(NA_real_, NA_real_))
  }

  one_pcor <- function(tested, other) {
    out <- try(
      psych::partial.r(data = d, x = c("x", tested), y = other,
                       use = pcor.use, method = cor.type)[2],
      silent = TRUE
    )
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
  }

  c(one_pcor("y", "z"), one_pcor("z", "y"))
}


#' Lower and upper limit of a boot.ci() result
#'
#' @param ci_int a result of boot::boot.ci() or a try-error
#' @param boot_ci_type the interval type that was requested
#'
#' @return a numeric vector of length 2, NA when no interval is available
#' @keywords internal
#' @noRd

boot_ci_bounds <- function(ci_int, boot_ci_type) {

  if (inherits(ci_int, "try-error") || is.null(ci_int)) {
    return(c(NA_real_, NA_real_))
  }

  bounds <- switch(boot_ci_type,
                   norm = ci_int$normal[2:3],
                   perc = ci_int$percent[4:5],
                   basic = ci_int$basic[4:5],
                   bca = ci_int$bca[4:5],
                   NULL)

  if (length(bounds) != 2) {
    return(c(NA_real_, NA_real_))
  }

  as.numeric(bounds)
}


#' Partial correlation of the response with the control window
#'
#' Internal helper of the seascorr functions with pcor_control = TRUE. For one
#' window, it calculates the partial correlation of the response with the
#' control variable, controlling for the primary variable, and its bootstrap
#' interval, in the same way as the functions calculate the primary partial
#' correlation. The result is what the functions would return for this window
#' with the primary and control data swapped. With boot() the second statistic
#' of the primary bootstrap run is used, so that no further random numbers
#' are drawn.
#'
#' @param data the window data, with columns x (response), y (primary) and z
#' (control), for a window with enough complete years
#' @param calc the boot() result of the primary calculation, computed with
#' boot_f_pcor_both(), or NULL
#'
#' @return a numeric vector: estimate, lower and upper limit
#' @keywords internal
#' @noRd

seascorr_control_stats <- function(data, boot, boot_fast_obj, calc,
                                   pcor_method, pcor_na_use, boot_conf_int,
                                   boot_ci_type) {

  out <- c(NA_real_, NA_real_, NA_real_)

  observed <- function() {
    as.numeric(psych::partial.r(data = data, x = c("x", "z"), y = "y",
                                use = pcor_na_use, method = pcor_method)[2])
  }

  if (!isTRUE(boot)) {
    out[1] <- observed()
    return(out)
  }

  if (!is.null(boot_fast_obj)) {

    boot_t <- boot_fast_pcor_replicates(boot_fast_obj, data[, "z"], data[, "y"])
    boot_t0 <- observed()

    return(c(boot_t0, boot_fast_ci(boot_t, boot_t0, conf = boot_conf_int,
                                   type = boot_ci_type)[1:2]))
  }

  ci_int <- try(boot::boot.ci(calc, conf = boot_conf_int, type = boot_ci_type,
                              index = 2L), silent = TRUE)

  c(as.numeric(calc$t0[2]), boot_ci_bounds(ci_int, boot_ci_type))
}


#' Significance filter of a seascorr surface
#'
#' Internal helper of daily_response_seascorr() and
#' monthly_response_seascorr(). It applies remove_insignificant_method to one
#' surface of partial correlations between the response and the windows in
#' climate_tested, each controlling for the paired window in climate_other.
#' It is called for the primary surface and, with pcor_control = TRUE, for the
#' control surface with the two lists swapped, so both surfaces are filtered
#' in exactly the same way. The permutation methods draw their permutations
#' once and return them, so that the second call can reuse them: both surfaces
#' are tested against the same permutations, and the control surface draws no
#' further random numbers.
#'
#' @param surface,surface_lower,surface_upper the surface and its bootstrap
#' limits
#' @param climate_tested,climate_other lists of aggregated windows, as
#' collected in the window loops
#' @param significance_filter the list describing the filter, to be completed
#' @param perm_index NULL, or permutations returned by an earlier call
#' @param label "primary" or "control", for the warnings
#'
#' @return a list with the filtered surface and limits, the completed
#' significance_filter and the permutations used (perm_index)
#' @keywords internal
#' @noRd

seascorr_significance <- function(surface, surface_lower, surface_upper,
                                  climate_tested, climate_other, response,
                                  significance_filter,
                                  remove_insignificant_method, alpha,
                                  pcor_method, pcor_na_use, perm_n,
                                  cluster_threshold, reference_window,
                                  perm_index = NULL, label = "primary") {

  # The permutations are drawn in the same order as sample(y_perm) in a loop,
  # so the results are those of earlier versions.
  draw_permutations <- function(n_years) {
    if (!is.null(perm_index) && nrow(perm_index) == n_years &&
        ncol(perm_index) == perm_n) {
      return(perm_index)
    }
    vapply(seq_len(perm_n), function(perm_i) sample.int(n_years),
           integer(n_years))
  }

  # Partial correlation controls for one variable, which the critical
  # value and p value helpers account for through df_adjust.
  n_obs <- nrow(response)
  original_non_missing <- !is.na(surface)
  keep <- NULL

  # 1 Plain alpha threshold
  if (remove_insignificant_method == "alpha") {

    critical_threshold_cor <- critical_r(n_obs, alpha = alpha,
                                         cor_method = pcor_method)

    keep <- abs(surface) >= abs(critical_threshold_cor)

    significance_filter$applied <- TRUE
    significance_filter$method <- "alpha"
    significance_filter$criterion <- paste0("|partial r| >= ",
                                            round(critical_threshold_cor, 3),
                                            " (alpha = ", alpha, ")")

    significance_filter$critical_r <- critical_threshold_cor

    # 2 Bootstrap confidence interval
  } else if (remove_insignificant_method == "boot") {

    if (is.null(surface_lower) ||
        is.null(surface_upper) ||
        !all(dim(surface) == dim(surface_lower)) ||
        !all(dim(surface) == dim(surface_upper))) {

      warning(paste0("remove_insignificant_method = 'boot' was requested, but ",
                     "bootstrap confidence-interval matrices are missing or ",
                     "do not match the calculation matrix. The significance ",
                     "filter is ignored."), call. = FALSE)

    } else {

      keep <- (surface_lower > 0 & surface_upper > 0) |
        (surface_lower < 0 & surface_upper < 0)

      significance_filter$applied <- TRUE
      significance_filter$method <- "boot"
      significance_filter$criterion <- "Bootstrap confidence interval excludes zero"
    }

    # 3 False discovery rate
  } else if (remove_insignificant_method == "fdr_BY" |
             remove_insignificant_method == "fdr_BH") {

    fdr_type <- sub("^fdr_", "", remove_insignificant_method)

    p_surface <- p_from_r(surface, n_obs,
                          cor_method = pcor_method, df_adjust = 1)

    p_adjusted <- stats::p.adjust(as.vector(p_surface), method = fdr_type)
    p_adjusted <- matrix(p_adjusted, nrow = nrow(p_surface))

    keep <- p_adjusted < alpha

    significance_filter$applied <- TRUE
    significance_filter$method <- remove_insignificant_method
    significance_filter$criterion <- paste0("Adjusted p < ", alpha,
                                            " (", fdr_type, ")")

    # 4 Max-statistic permutation
  } else if (remove_insignificant_method == "maxstat") {

    climate_x1 <- do.call(cbind, climate_tested)
    climate_x2 <- do.call(cbind, climate_other)

    if (is.null(climate_x1) || is.null(climate_x2) ||
        ncol(climate_x1) < 1 ||
        ncol(climate_x1) != ncol(climate_x2)) {

      warning(paste0("remove_insignificant_method = 'maxstat' was requested, ",
                     "but the matrices of aggregated primary and control ",
                     "windows are not available or do not match. The ",
                     "significance filter is ignored."), call. = FALSE)

    } else {

      shared_years <- intersect(row.names(response), row.names(climate_x1))

      if (length(shared_years) < 10) {

        warning(paste0("remove_insignificant_method = 'maxstat' was requested, ",
                       "but fewer than 10 years are shared between the response ",
                       "and the aggregated climate windows. The significance ",
                       "filter is ignored."), call. = FALSE)

      } else {

        y_perm  <- response[shared_years, 1]
        x1_perm <- climate_x1[shared_years, , drop = FALSE]
        x2_perm <- climate_x2[shared_years, , drop = FALSE]

        # Correlation between each tested window and its paired window of the
        # other variable. This does not involve the response, so it is
        # unchanged by permutation and is computed only once.
        r_yz <- vapply(seq_len(ncol(x1_perm)), function(win_i) {
          suppressWarnings(stats::cor(x1_perm[, win_i], x2_perm[, win_i],
                                      method = pcor_method, use = pcor_na_use))
        }, numeric(1))

        denom_z <- sqrt(1 - r_yz ^ 2)

        perm_index <- draw_permutations(length(y_perm))

        # Null distribution of the largest absolute partial correlation
        # obtained anywhere in the search when the response is unrelated to
        # climate. The same shuffle is applied to both correlations so that
        # the response stays internally consistent within a permutation.
        null_max <- rep(NA_real_, perm_n)

        for (perm_i in seq_len(perm_n)) {

          y_shuffled <- y_perm[perm_index[, perm_i]]

          r_xy <- suppressWarnings(
            as.numeric(stats::cor(y_shuffled, x1_perm,
                                  method = pcor_method, use = pcor_na_use)))

          r_xz <- suppressWarnings(
            as.numeric(stats::cor(y_shuffled, x2_perm,
                                  method = pcor_method, use = pcor_na_use)))

          perm_pcor <- (r_xy - r_xz * r_yz) / (sqrt(1 - r_xz ^ 2) * denom_z)

          null_max[perm_i] <- max(abs(perm_pcor), na.rm = TRUE)
        }

        critical_threshold_perm <- as.numeric(
          stats::quantile(null_max, probs = 1 - alpha, na.rm = TRUE))

        observed_max <- max(abs(surface), na.rm = TRUE)

        keep <- abs(surface) >= critical_threshold_perm

        p_global <- (sum(null_max >= observed_max, na.rm = TRUE) + 1) /
          (perm_n + 1)

        significance_filter$applied <- TRUE
        significance_filter$method <- "maxstat"
        significance_filter$criterion <- paste0("|partial r| >= ",
                                                round(critical_threshold_perm, 3),
                                                " (max-statistic permutation, ",
                                                perm_n, " permutations, alpha = ",
                                                alpha, ")")
        significance_filter$critical_r <- critical_threshold_perm
        significance_filter$observed_max_r <- observed_max
        significance_filter$p_global <- p_global
        significance_filter$perm_n <- perm_n
        significance_filter$n_windows_searched <- ncol(x1_perm)
      }
    }

    # 5 Sidak correction on the effective number of tests
  } else if (remove_insignificant_method == "sidak") {

    climate_x1 <- do.call(cbind, climate_tested)
    climate_x2 <- do.call(cbind, climate_other)

    if (is.null(climate_x1) || is.null(climate_x2) ||
        ncol(climate_x1) < 2 ||
        ncol(climate_x1) != ncol(climate_x2)) {

      warning(paste0("remove_insignificant_method = 'sidak' was requested, ",
                     "but the matrices of aggregated primary and control ",
                     "windows are not available or do not match. The ",
                     "significance filter is ignored."), call. = FALSE)

    } else {

      shared_years <- intersect(row.names(response), row.names(climate_x1))

      if (length(shared_years) < 10) {

        warning(paste0("remove_insignificant_method = 'sidak' was requested, ",
                       "but fewer than 10 years are shared between the ",
                       "response and the aggregated climate windows. The ",
                       "significance filter is ignored."), call. = FALSE)

      } else {

        n_obs_sidak <- length(shared_years)

        x1_sidak <- climate_x1[shared_years, , drop = FALSE]
        x2_sidak <- climate_x2[shared_years, , drop = FALSE]

        # Missing years are left out rather than whole windows, so that the
        # effective number of tests describes every screened window, see
        # sidak_window_matrix(). A year or window counts as missing when the
        # value of either variable is missing. Complete data are not
        # affected.
        x1_sidak <- sidak_window_matrix(x1_sidak, x_control = x2_sidak)
        rm(x2_sidak)

        if (ncol(x1_sidak) > 1) {
          x1_sidak <- x1_sidak[, apply(x1_sidak, 2, stats::sd) > 0,
                               drop = FALSE]
        }

        if (ncol(x1_sidak) < 2) {

          warning(paste0("remove_insignificant_method = 'sidak' was ",
                         "requested, but fewer than two ", label, " windows ",
                         "vary. The significance filter is ignored."),
                  call. = FALSE)

        } else {

          # Effective number of independent tests among the screened windows,
          # following Li and Ji (2005). The correlation matrix among windows
          # has rank at most the number of years, so its non-zero eigenvalues
          # come from the years-by-years matrix rather than the
          # windows-by-windows one.
          z_sidak <- scale(x1_sidak)
          gram <- tcrossprod(z_sidak) / (nrow(z_sidak) - 1)

          eig <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
          eig <- eig[eig > 1e-10]

          m_eff <- sum((eig >= 1) + (eig - floor(eig)))

          # One degree of freedom is spent on the controlled variable
          alpha_sidak <- 1 - (1 - alpha) ^ (1 / m_eff)
          critical_threshold_sidak <- critical_r(n_obs_sidak - 1,
                                                 alpha = alpha_sidak,
                                                 cor_method = pcor_method)

          keep <- abs(surface) >= critical_threshold_sidak

          significance_filter$applied <- TRUE
          significance_filter$method <- "sidak"
          significance_filter$criterion <- paste0("|partial r| >= ",
                                                  round(critical_threshold_sidak, 3),
                                                  " (Sidak on ", round(m_eff, 1),
                                                  " effective tests, alpha = ",
                                                  alpha, ")")
          significance_filter$critical_r <- critical_threshold_sidak
          significance_filter$effective_tests <- m_eff
          significance_filter$alpha_per_test <- alpha_sidak
          significance_filter$n_windows_searched <- ncol(x1_sidak)
        }
      }
    }

    # 6 Cluster-extent inference
  } else if (remove_insignificant_method == "cluster") {

    climate_x1 <- do.call(cbind, climate_tested)
    climate_x2 <- do.call(cbind, climate_other)

    if (is.null(climate_x1) || is.null(climate_x2) ||
        ncol(climate_x1) < 1 ||
        ncol(climate_x1) != ncol(climate_x2)) {

      warning(paste0("remove_insignificant_method = 'cluster' was requested, ",
                     "but the matrices of aggregated primary and control ",
                     "windows are not available or do not match. The ",
                     "significance filter is ignored."), call. = FALSE)

    } else {

      shared_years <- intersect(row.names(response), row.names(climate_x1))

      if (length(shared_years) < 10) {

        warning(paste0("remove_insignificant_method = 'cluster' was ",
                       "requested, but fewer than 10 years are shared ",
                       "between the response and the aggregated climate ",
                       "windows. The significance filter is ignored."),
                call. = FALSE)

      } else {

        y_perm <- response[shared_years, 1]
        x1_perm <- climate_x1[shared_years, , drop = FALSE]
        x2_perm <- climate_x2[shared_years, , drop = FALSE]
        n_obs_cl <- length(shared_years)

        # Cluster inference needs the permuted statistics laid out on the
        # same surface as the observed ones, unlike the max-statistic test
        # which needs only their maximum. Each aggregated window is named
        # "start_end", which gives its width and position; the mapping is
        # checked against the observed surface before it is used.
        win_nm <- colnames(x1_perm)
        parts <- strsplit(win_nm, "_", fixed = TRUE)
        w_start <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 1)))
        w_end <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
        w_width <- w_end - w_start + 1

        surf_widths <- suppressWarnings(as.numeric(row.names(surface)))

        if (nrow(surface) == 1) {
          row_idx <- rep(1L, length(w_width))
        } else {
          row_idx <- match(w_width, surf_widths)
        }

        col_idx <- if (reference_window == "start") {
          w_start
        } else if (reference_window == "end") {
          w_end
        } else {
          round2(w_start + w_width / 2, 0)
        }

        usable <- !is.na(row_idx) & !is.na(col_idx) &
          row_idx >= 1 & row_idx <= nrow(surface) &
          col_idx >= 1 & col_idx <= ncol(surface)

        map_ok <- FALSE
        if (any(usable)) {

          r_xy_o <- suppressWarnings(
            as.numeric(stats::cor(y_perm, x1_perm[, usable, drop = FALSE],
                                  method = pcor_method, use = pcor_na_use)))
          r_xz_o <- suppressWarnings(
            as.numeric(stats::cor(y_perm, x2_perm[, usable, drop = FALSE],
                                  method = pcor_method, use = pcor_na_use)))
          r_yz_o <- vapply(which(usable), function(win_i) {
            suppressWarnings(stats::cor(x1_perm[, win_i], x2_perm[, win_i],
                                        method = pcor_method,
                                        use = pcor_na_use))
          }, numeric(1))

          obs_pcor <- (r_xy_o - r_xz_o * r_yz_o) /
            (sqrt(1 - r_xz_o ^ 2) * sqrt(1 - r_yz_o ^ 2))

          on_surface <- surface[cbind(row_idx[usable], col_idx[usable])]
          both <- is.finite(obs_pcor) & is.finite(on_surface)
          map_ok <- sum(both) > 0 &&
            max(abs(obs_pcor[both] - on_surface[both])) < 1e-8
        }

        if (!map_ok) {

          warning(paste0("remove_insignificant_method = 'cluster' was ",
                         "requested, but the aggregated windows could not be ",
                         "matched to their positions on the response ",
                         "surface. The significance filter is ignored."),
                  call. = FALSE)

        } else {

          idx_cells <- cbind(row_idx[usable], col_idx[usable])
          x1_use <- x1_perm[, usable, drop = FALSE]
          x2_use <- x2_perm[, usable, drop = FALSE]

          form_thr <- if (identical(cluster_threshold, "nominal")) {
            critical_r(n_obs_cl, alpha = alpha, cor_method = pcor_method)
          } else {
            cluster_threshold
          }

          # the correlation between the two variables does not involve the
          # response, so it is unchanged by permutation and computed only once
          r_yz <- vapply(seq_len(ncol(x1_use)), function(win_i) {
            suppressWarnings(stats::cor(x1_use[, win_i], x2_use[, win_i],
                                        method = pcor_method,
                                        use = pcor_na_use))
          }, numeric(1))
          denom_z <- sqrt(1 - r_yz ^ 2)

          perm_index <- draw_permutations(length(y_perm))

          null_clust <- rep(NA_real_, perm_n)
          null_max <- rep(NA_real_, perm_n)

          perm_surface <- matrix(NA_real_, nrow(surface), ncol(surface))

          for (perm_i in seq_len(perm_n)) {

            y_shuffled <- y_perm[perm_index[, perm_i]]

            r_xy <- suppressWarnings(
              as.numeric(stats::cor(y_shuffled, x1_use,
                                    method = pcor_method, use = pcor_na_use)))
            r_xz <- suppressWarnings(
              as.numeric(stats::cor(y_shuffled, x2_use,
                                    method = pcor_method, use = pcor_na_use)))

            perm_pcor <- (r_xy - r_xz * r_yz) /
              (sqrt(1 - r_xz ^ 2) * denom_z)

            perm_surface[] <- NA_real_
            perm_surface[idx_cells] <- perm_pcor

            null_max[perm_i] <- max(abs(perm_pcor), na.rm = TRUE)
            null_clust[perm_i] <- max(cluster_fast_masses(perm_surface,
                                                          form_thr))
          }

          crit_mass <- as.numeric(
            stats::quantile(null_clust, probs = 1 - alpha, na.rm = TRUE))
          crit_max_cl <- as.numeric(
            stats::quantile(null_max, probs = 1 - alpha, na.rm = TRUE))

          obs_mask <- !is.na(surface) & abs(surface) > form_thr
          cc <- cluster_fast_label(obs_mask)

          keep <- matrix(FALSE, nrow(surface), ncol(surface))
          n_retained <- 0L
          obs_masses <- numeric(0)

          if (cc$n > 0L) {
            for (cl in seq_len(cc$n)) {
              in_cl <- cc$labels == cl
              mass <- sum(abs(surface[in_cl]), na.rm = TRUE)
              obs_masses <- c(obs_masses, mass)
              if (mass > crit_mass) {
                keep[in_cl] <- TRUE
                n_retained <- n_retained + 1L
              }
            }
          }

          significance_filter$applied <- TRUE
          significance_filter$method <- "cluster"
          significance_filter$criterion <- paste0("cluster mass > ",
                                                  round(crit_mass, 1),
                                                  ", clusters formed at ",
                                                  "|partial r| > ",
                                                  round(form_thr, 3), " (",
                                                  perm_n,
                                                  " permutations, alpha = ",
                                                  alpha, ")")
          significance_filter$forming_threshold <- form_thr
          significance_filter$cluster_mass_critical <- crit_mass
          significance_filter$n_clusters_formed <- cc$n
          significance_filter$n_clusters_retained <- n_retained
          significance_filter$largest_cluster_mass <-
            if (length(obs_masses)) max(obs_masses) else 0
          significance_filter$maxstat_critical_r <- crit_max_cl
          significance_filter$perm_n <- perm_n
          significance_filter$n_windows_searched <- ncol(x1_use)
        }
      }
    }
  }

  # Apply whichever mask was produced
  if (!is.null(keep)) {

    keep[is.na(keep)] <- FALSE

    surface[!keep] <- NA

    # The same mask is applied to the confidence-interval matrices, so that
    # a filtered surface is not accompanied by intervals for windows the
    # filter removed.
    if (!is.null(surface_lower) &&
        all(dim(surface_lower) == dim(surface))) {
      surface_lower[!keep] <- NA
      surface_upper[!keep] <- NA
    }

    significance_filter$removed_cells <- sum(original_non_missing & !keep,
                                             na.rm = TRUE)
    significance_filter$retained_cells <- sum(!is.na(surface))
    significance_filter$total_cells <- length(surface)
  }

  list(surface = surface,
       lower = surface_lower,
       upper = surface_upper,
       significance_filter = significance_filter,
       perm_index = perm_index)
}


#' Both partial correlations and the shared variance of the optimal window
#'
#' Internal helper of the seascorr functions. From the response, the primary
#' and the control series of the optimal window, it returns the three simple
#' correlations, the two partial correlations and how the variance of the
#' response explained by both variables together (R squared of a regression on
#' both) splits into a part unique to the primary variable, a part unique to
#' the control variable and a part they share. The unique parts are squared
#' semipartial correlations. The shared part is negative when one variable
#' suppresses the other. The split is based on the correlations of pcor_method
#' and pcor_na_use, so for "spearman" it refers to ranks. For "kendall", which
#' does not decompose variance, it is NA.
#'
#' @param data a matrix or data frame with the response, the primary and the
#' control series in this order
#'
#' @return a data frame with the columns Variable and Value
#' @keywords internal
#' @noRd

seascorr_variance_partitioning <- function(data, pcor_method, pcor_na_use) {

  data <- as.matrix(data[, 1:3, drop = FALSE])
  colnames(data) <- c("x", "y", "z")

  r <- suppressWarnings(stats::cor(data, method = pcor_method,
                                   use = pcor_na_use))
  r_xy <- r[1, 2]
  r_xz <- r[1, 3]
  r_yz <- r[2, 3]

  # the same calculation as for the surfaces, so the values agree exactly
  partial <- function(tested, other) {
    out <- try(as.numeric(psych::partial.r(data = data, x = c("x", tested),
                                           y = other, use = pcor_na_use,
                                           method = pcor_method)[2]),
               silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else out
  }

  if (identical(pcor_method, "kendall")) {

    r_squared <- unique_primary <- unique_control <- shared <- NA_real_

  } else {

    r_squared <- (r_xy ^ 2 + r_xz ^ 2 - 2 * r_xy * r_xz * r_yz) / (1 - r_yz ^ 2)
    unique_primary <- r_squared - r_xz ^ 2
    unique_control <- r_squared - r_xy ^ 2
    shared <- r_squared - unique_primary - unique_control
  }

  data.frame(Variable = c("r_response_primary",
                          "r_response_control",
                          "r_primary_control",
                          "pcor_primary",
                          "pcor_control",
                          "r_squared",
                          "unique_primary",
                          "unique_control",
                          "shared"),
             Value = as.numeric(c(r_xy, r_xz, r_yz,
                                  partial("y", "z"), partial("z", "y"),
                                  r_squared, unique_primary, unique_control,
                                  shared)),
             stringsAsFactors = FALSE)
}
