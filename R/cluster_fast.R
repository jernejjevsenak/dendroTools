#' cluster_fast_label
#'
#' Label connected components of a logical matrix using four-connectivity.
#'
#' Implemented by run-length encoding each row and then joining runs that
#' overlap in adjacent rows, rather than visiting every cell. The number of
#' runs is typically orders of magnitude smaller than the number of cells,
#' which matters because cluster-extent inference has to repeat this for every
#' permutation.
#'
#' @param mask logical matrix
#'
#' @return a list with the integer label matrix and the number of clusters
#' @keywords internal
#'
cluster_fast_label <- function(mask) {

  nr <- nrow(mask)
  nc <- ncol(mask)

  runs <- vector("list", nr)
  n_lab <- 0L

  for (i in seq_len(nr)) {

    v <- mask[i, ]
    v[is.na(v)] <- FALSE
    if (!any(v)) next

    r <- rle(v)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    sel <- r$values

    s <- starts[sel]
    e <- ends[sel]
    lab <- n_lab + seq_along(s)
    n_lab <- n_lab + length(s)

    runs[[i]] <- cbind(start = s, end = e, lab = lab, row = i)
  }

  if (n_lab == 0L) {
    return(list(labels = matrix(0L, nr, nc), n = 0L))
  }

  # union-find over run labels
  parent <- seq_len(n_lab)

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }

  join <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[ra] <<- rb
    invisible(NULL)
  }

  # runs in consecutive rows belong to the same cluster when their column
  # ranges overlap
  for (i in seq_len(nr - 1L)) {

    A <- runs[[i]]
    B <- runs[[i + 1L]]
    if (is.null(A) || is.null(B)) next

    for (a in seq_len(nrow(A))) {
      lo_a <- A[a, "start"]; hi_a <- A[a, "end"]
      hit <- which(B[, "start"] <= hi_a & B[, "end"] >= lo_a)
      for (b in hit) join(A[a, "lab"], B[b, "lab"])
    }
  }

  # resolve to consecutive cluster numbers
  roots <- vapply(seq_len(n_lab), find_root, integer(1))
  uniq <- sort(unique(roots))
  remap <- integer(max(uniq))
  remap[uniq] <- seq_along(uniq)

  labels <- matrix(0L, nr, nc)
  for (i in seq_len(nr)) {
    A <- runs[[i]]
    if (is.null(A)) next
    for (a in seq_len(nrow(A))) {
      labels[i, A[a, "start"]:A[a, "end"]] <- remap[roots[A[a, "lab"]]]
    }
  }

  list(labels = labels, n = length(uniq))
}


#' cluster_fast_masses
#'
#' Mass of every supra-threshold cluster on a statistic surface. Mass is the
#' sum of absolute values within a cluster, so it responds to both the size of
#' a region and the strength of the values inside it.
#'
#' @param r_mat numeric matrix of the statistic
#' @param thr forming threshold on the absolute statistic
#'
#' @return numeric vector of cluster masses, or 0 when nothing forms
#' @keywords internal
#'
cluster_fast_masses <- function(r_mat, thr) {

  mask <- !is.na(r_mat) & abs(r_mat) > thr
  if (!any(mask)) return(0)

  cc <- cluster_fast_label(mask)
  if (cc$n == 0L) return(0)

  a <- abs(r_mat)
  a[is.na(a)] <- 0

  as.numeric(tapply(a[cc$labels > 0], cc$labels[cc$labels > 0], sum))
}
