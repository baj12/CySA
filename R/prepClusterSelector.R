#' Prepare Data for the Cluster Selector Shiny App
#'
#' Subsamples a \code{SingleCellExperiment} object and builds the inputs
#' required by \code{\link{clusterSelector}}.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing the full dataset.
#' @param somFile Path to the SOM object. Used to cache/restore the subsampled
#'   version.
#' @param dList Optional list of marker pairs for 2D plots. If \code{NULL},
#'   defaults to the first 12 row names of \code{sce}.
#' @param total_cells_to_sample Number of cells to subsample in total.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with \code{sce}, \code{sce_subsampled}, and \code{dList}.
#'
#' @examples
#' sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
#' prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
#' names(prepped)
#'
#' @export
prepClusterSelectorData <- function(sce,
                                    somFile = NULL,
                                    dList = NULL,
                                    total_cells_to_sample = 100000,
                                    somCodesName = "SOM_codes",
                                    assay = "exprs",
                                    seed = 123) {

  if (is.null(dList)) {
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed
    rn <- if (!is.null(colsUsed) && length(colsUsed) >= 12) colsUsed else rownames(sce)
    if (length(rn) < 12) {
      stop("sce must have at least 12 row names to build default dList")
    }
    dList <- list(
      d1 = c(rn[1], rn[2]),
      d2 = c(rn[3], rn[4]),
      d3 = c(rn[5], rn[6]),
      d4 = c(rn[7], rn[8]),
      d5 = c(rn[9], rn[10]),
      d6 = c(rn[11], rn[12])
    )
  }

  # Cache file handling
  cache_file <- NULL
  if (!is.null(somFile)) {
    cache_file <- paste0(tools::file_path_sans_ext(somFile), ".subsampled.RData")
  }

  if (!is.null(cache_file) && file.exists(cache_file)) {
    env <- new.env()
    load(cache_file, envir = env)
    sce_subsampled <- env$sce_subsampled
  } else {
    set.seed(seed)

    cd <- as.data.frame(SingleCellExperiment::colData(sce))
    proportions_df <- dplyr::group_by(cd, .data$sample_id, .data$cluster_id)
    proportions_df <- dplyr::summarise(proportions_df, group_size = dplyr::n(), .groups = "drop")
    proportions_df <- dplyr::ungroup(proportions_df)
    proportions_df <- dplyr::mutate(
      proportions_df,
      total_size = sum(.data$group_size),
      proportion = .data$group_size / .data$total_size,
      n_to_sample = ceiling(.data$proportion * total_cells_to_sample + 10)
    )

    sampling_indices <- purrr::map2(
      proportions_df$n_to_sample,
      proportions_df$group_size,
      ~{
        if (.x > .y) return(seq(.y))
        sample(.y, .x)
      }
    )
    sampling_indices <- purrr::set_names(
      sampling_indices,
      paste(proportions_df$sample_id, proportions_df$cluster_id, sep = "_")
    )
    sampling_indices <- purrr::imap(sampling_indices, ~{
      sid <- stringr::str_replace(.y, "(^.*)_(.*)", "\\1")
      cid <- stringr::str_replace(.y, "(^.*)_(.*)", "\\2")
      sce[,
        SingleCellExperiment::colData(sce)$sample_id == sid &
          SingleCellExperiment::colData(sce)$cluster_id == cid
      ][, .x]
    })

    sce_subsampled <- do.call(SummarizedExperiment::cbind, sampling_indices)

    if (!is.null(cache_file)) {
      save(sce_subsampled, file = cache_file)
    }
  }

  if (length(SummarizedExperiment::assays(sce_subsampled)) == 0) {
    SummarizedExperiment::assays(sce_subsampled)[[1]] <-
      SummarizedExperiment::assays(sce)[[1]][, SummarizedExperiment::colData(sce_subsampled)$id]
    names(SummarizedExperiment::assays(sce_subsampled)) <- assay
  }

  list(sce = sce, sce_subsampled = sce_subsampled, dList = dList)
}
