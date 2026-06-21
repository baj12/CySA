# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

#' Minimal example SingleCellExperiment for CySA
#'
#' Creates a small, deterministic \code{SingleCellExperiment} object that
#' contains the metadata and column data expected by \code{CySA} functions.
#'
#' @param n_cells Number of cells (columns).
#' @param n_nodes Number of SOM nodes.
#' @param n_markers Number of markers (rows).
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with
#'   \code{SOM_codes}, \code{SOM_stats}, \code{sample_id}, and
#'   \code{cluster_id}.
#'
#' @examples
#' sce <- CySA_example_sce()
#' head(S4Vectors::metadata(sce)$SOM_codes)
#'
#' @export
CySA_example_sce <- function(n_cells = 1000, n_nodes = 50, n_markers = 12) {
  markers <- paste0("marker", seq_len(n_markers))

  counts <- matrix(
    data = seq_len(n_cells * n_markers) %% 1000L + 1L,
    nrow = n_markers,
    ncol = n_cells
  )
  rownames(counts) <- markers
  colnames(counts) <- paste0("cell", seq_len(n_cells))

  sample_ids <- rep(paste0("sample", seq_len(4)), length.out = n_cells)
  cluster_ids <- rep(seq_len(n_nodes), length.out = n_cells)
  cd <- S4Vectors::DataFrame(
    sample_id = sample_ids,
    cluster_id = cluster_ids
  )

  som_codes <- matrix(
    data = seq_len(n_nodes * n_markers) / (n_nodes * n_markers),
    nrow = n_nodes,
    ncol = n_markers
  )
  colnames(som_codes) <- markers
  rownames(som_codes) <- as.character(seq_len(n_nodes))

  som_stats <- data.frame(
    id = seq_len(n_nodes),
    n = as.integer(table(factor(cluster_ids, levels = seq_len(n_nodes)))),
    mean = seq_len(n_nodes) / n_nodes,
    median = seq_len(n_nodes) / n_nodes,
    rdQu = seq_len(n_nodes) / n_nodes,
    max = seq_len(n_nodes) / n_nodes,
    stringsAsFactors = FALSE
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(exprs = counts),
    colData = cd
  )

  S4Vectors::metadata(sce)$SOM_codes <- som_codes
  S4Vectors::metadata(sce)$SOM_stats <- som_stats
  S4Vectors::metadata(sce)$map <- list(colsUsed = markers)

  experiment_info <- data.frame(
    sample_id = unique(sample_ids),
    some_numeric = seq_along(unique(sample_ids)),
    stringsAsFactors = FALSE
  )
  S4Vectors::metadata(sce)$experiment_info <- experiment_info

  SingleCellExperiment::int_metadata(sce)$channels <- S4Vectors::DataFrame(
    channel_name = markers,
    marker_name = markers
  )

  sce
}
