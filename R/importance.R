# ==========================================================================
# Reliability importance measures
# ==========================================================================
#
# Three classical measures, all derivable from phi, min_cuts, and the
# component distributions:
#
#   - Birnbaum (reliability) importance:   dR/dp_j
#   - Criticality (Fussell):               component j failed AND critical | system failed
#   - Vesely-Fussell:                      some cut containing j is fully failed | system failed
#
# All live as default methods on `dist_structure`. Specialized subclasses
# may override for closed-form speed.
# ==========================================================================


#' @rdname birnbaum_importance
#' @export
birnbaum_importance.dist_structure <- function(x, j, p) {
  m <- ncomponents(x)
  if (length(p) == 1L) p <- rep(p, m)
  stopifnot(length(p) == m, all(p >= 0 & p <= 1),
            j >= 1L, j <= m, length(j) == 1L)
  p1 <- p; p1[j] <- 1
  p0 <- p; p0[j] <- 0
  reliability(x, p1) - reliability(x, p0)
}


# Helper: component survivals at time t, as a length-m numeric vector.
component_survivals_at <- function(x, t) {
  m <- ncomponents(x)
  vapply(seq_len(m), function(k) {
    algebraic.dist::surv(component(x, k))(t)
  }, numeric(1L))
}


#' @rdname criticality_importance
#' @export
criticality_importance.dist_structure <- function(x, j, t) {
  m <- ncomponents(x)
  stopifnot(j >= 1L, j <= m, length(j) == 1L)
  S <- component_survivals_at(x, t)
  F_j <- 1 - S[j]
  F_sys <- 1 - reliability(x, S)
  if (F_sys == 0) return(0)
  ib <- birnbaum_importance(x, j, S)
  as.numeric(ib * F_j / F_sys)
}


#' @rdname vesely_fussell_importance
#' @export
vesely_fussell_importance.dist_structure <- function(x, j, t) {
  m <- ncomponents(x)
  stopifnot(j >= 1L, j <= m, length(j) == 1L)
  cuts <- min_cuts(x)
  cuts_j <- Filter(function(K) j %in% K, cuts)
  if (length(cuts_j) == 0L) return(0)
  S <- component_survivals_at(x, t)
  F_vec <- 1 - S
  F_sys <- 1 - reliability(x, S)
  if (F_sys == 0) return(0)
  # Inclusion-exclusion: P(union of events "cut K_l is fully failed")
  # where each event's probability under independence is prod_{i in K_l} F_i.
  # For an intersection of several cut-failure events, all components in
  # the union of those cuts must have failed:
  # P(cut A failed AND cut B failed) = prod_{i in A union B} F_i.
  nC <- length(cuts_j)
  p_union <- 0
  for (sz in seq_len(nC)) {
    subsets <- utils::combn(nC, sz, simplify = FALSE)
    for (subset in subsets) {
      union_idx <- sort(unique(unlist(cuts_j[subset])))
      p_union <- p_union + (-1)^(sz + 1L) * prod(F_vec[union_idx])
    }
  }
  as.numeric(p_union / F_sys)
}
