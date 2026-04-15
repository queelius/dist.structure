# ==========================================================================
# Compositional operations
# ==========================================================================
#
# Build new dist_structure objects from existing ones without leaving the
# protocol:
#
#   - substitute_component(x, j, new_component): component-level edit
#   - compose_systems(outer, inner_list):        hierarchical nesting
# ==========================================================================


#' @rdname substitute_component
#' @export
substitute_component.dist_structure <- function(x, j, new_component) {
  m <- ncomponents(x)
  stopifnot(j >= 1L, j <= m, length(j) == 1L)
  # Collect existing components via the primitive, so this works for any
  # subclass providing component().
  components <- lapply(seq_len(m), function(k) component(x, k))
  components[[j]] <- new_component
  coherent_dist(
    min_paths = min_paths(x),
    components = components,
    m = m
  )
}


# Helper: return the component list of an inner argument (a dist_structure
# or a plain dist), along with its min_paths treated relative to a
# zero-based offset of its first component.
unpack_inner <- function(inner) {
  if (is_dist_structure(inner)) {
    n <- ncomponents(inner)
    components <- lapply(seq_len(n), function(k) component(inner, k))
    paths <- min_paths(inner)
  } else {
    # Plain dist: treated as a single-component series (its only path is
    # itself).
    components <- list(inner)
    paths <- list(1L)
  }
  list(components = components, paths = paths, size = length(components))
}


#' @rdname compose_systems
#' @export
compose_systems.dist_structure <- function(outer, inner_list) {
  m_outer <- ncomponents(outer)
  stopifnot(length(inner_list) == m_outer)
  unpacked <- lapply(inner_list, unpack_inner)
  sizes <- vapply(unpacked, function(u) u$size, integer(1L))
  offsets <- c(0L, cumsum(sizes))[seq_len(m_outer)]
  # Flatten components in outer order.
  all_components <- do.call(c, lapply(unpacked, function(u) u$components))
  # For each inner, shift its min_paths by its offset to get global indices.
  shifted_paths <- lapply(seq_len(m_outer), function(k) {
    lapply(unpacked[[k]]$paths, function(P) as.integer(P + offsets[k]))
  })
  # Composed min_paths: for each outer path P_out, take the Cartesian
  # product of inner paths for k in P_out and union their (shifted) indices.
  composed_paths <- list()
  for (P_out in min_paths(outer)) {
    inner_choices <- shifted_paths[P_out]
    grid <- expand.grid(lapply(inner_choices, seq_along),
                        KEEP.OUT.ATTRS = FALSE)
    for (row_idx in seq_len(nrow(grid))) {
      selected <- Map(`[[`, inner_choices, as.integer(grid[row_idx, ]))
      composed_paths[[length(composed_paths) + 1L]] <-
        sort(unique(unlist(selected)))
    }
  }
  composed_paths <- minimize_sets(composed_paths)
  coherent_dist(
    min_paths = composed_paths,
    components = all_components,
    m = sum(sizes)
  )
}
