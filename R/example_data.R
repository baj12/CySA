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
        cluster_id = factor(cluster_ids, levels = seq_len(n_nodes))
    )

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(exprs = counts),
        colData = cd
    )

    S4Vectors::metadata(sce)$SOM_codes <- .exampleSomCodes(n_nodes, markers)
    S4Vectors::metadata(sce)$SOM_stats <- .exampleSomStats(n_nodes,
        cluster_ids)
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

#' Deterministic SOM Codes of the Example Data
#'
#' @param n_nodes Number of SOM nodes.
#' @param markers Marker names.
#'
#' @return A nodes x markers matrix with increasing values in (0, 1].
#'
#' @keywords internal
.exampleSomCodes <- function(n_nodes, markers) {
    n <- n_nodes * length(markers)
    codes <- matrix(seq_len(n) / n, nrow = n_nodes, ncol = length(markers))
    colnames(codes) <- markers
    rownames(codes) <- as.character(seq_len(n_nodes))
    codes
}

#' Per-Node Statistics of the Example Data
#'
#' @param n_nodes Number of SOM nodes.
#' @param cluster_ids Node id of every cell.
#'
#' @return A data frame with \code{id}, cell count \code{n} and the
#'   \code{mean}, \code{median}, \code{rdQu}, \code{max} statistics.
#'
#' @keywords internal
.exampleSomStats <- function(n_nodes, cluster_ids) {
    ramp <- seq_len(n_nodes) / n_nodes
    data.frame(
        id = seq_len(n_nodes),
        n = as.integer(table(factor(cluster_ids, levels = seq_len(n_nodes)))),
        mean = ramp, median = ramp, rdQu = ramp, max = ramp,
        stringsAsFactors = FALSE
    )
}
