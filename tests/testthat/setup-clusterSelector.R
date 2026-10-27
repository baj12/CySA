# Setup helpers for clusterSelector tests.
#
# This file creates a reusable small example app and shared inputs. It is sourced
# once by testthat before any test file runs, so objects created here are
# available in all tests (but not exported from the package).

library(shiny)

# Build a small, deterministic app object that tests can share.
make_test_app <- function(n_cells = 200, n_nodes = 10, fsom = NULL) {
  set.seed(12)
  sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes)

  # Add a factor grouping column so the groupsVar/group1/group2 observers
  # (observers_clusterSelector.R:162-199) have something real to filter on.
  # experiment_info's built-in columns (sample_id, some_numeric) are
  # character/integer, never factors, so levels() on them is always NULL
  # and those observers' guard would otherwise never pass.
  S4Vectors::metadata(sce)$experiment_info$treatment <- factor(
    rep(c("control", "treated"), length.out = nrow(S4Vectors::metadata(sce)$experiment_info)),
    levels = c("control", "treated")
  )

  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
  som_codes <- S4Vectors::metadata(sce)$SOM_codes
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)),
    label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id
  )
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  }
  arr <- array(
    data = seq_len(10 * 10 * length(markers)),
    dim = c(10, 10, length(markers))
  )
  somRasterObj <- raster::brick(arr)
  names(somRasterObj) <- markers
  env <- new.env()
  app <- clusterSelector(
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList,
    dend = dend,
    dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = somRasterObj,
    fsom = fsom,
    env = env
  )
  list(
    app = app,
    env = env,
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled
  )
}

fsom_stub <- structure(
  list(MST = list(l = matrix(c(0, 1, 0, 1), ncol = 2))),
  class = "FlowSOM"
)

#' Test app with known cell-type labels (node_labels feature).
#'
#' Same fixture as make_test_app(), with per-node cell types and
#' cellTypeCol passed to prepClusterSelectorData so metadata(sce)$node_labels
#' exists.
make_labeled_test_app <- function(n_cells = 200, n_nodes = 10,
    cell_types = paste0("type", rep(1:5, length.out = n_nodes))) {
  set.seed(12)
  sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes,
    cell_types = cell_types)
  S4Vectors::metadata(sce)$experiment_info$treatment <- factor(
    rep(c("control", "treated"), length.out = nrow(S4Vectors::metadata(sce)$experiment_info)),
    levels = c("control", "treated")
  )
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  som_codes <- S4Vectors::metadata(sce)$SOM_codes
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)),
    label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id
  )
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  }
  arr <- array(
    data = seq_len(10 * 10 * length(markers)),
    dim = c(10, 10, length(markers))
  )
  somRasterObj <- raster::brick(arr)
  names(somRasterObj) <- markers
  env <- new.env()
  app <- clusterSelector(
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList,
    dend = dend,
    dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = somRasterObj,
    fsom = NULL,
    env = env
  )
  list(
    app = app,
    env = env,
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled
  )
}

# Two markers on different scales, a selection of nodes 1 and 2; marker
# "b" is listed first in the assay but separates less, so sorting reorders
# the rows relative to the matrix. Shared by the ridge test files.
make_interactive_sce <- function() {
  set.seed(3)
  n <- 200
  cl <- rep(c(1L, 2L, 3L, 4L), each = n / 4)
  a <- rnorm(n, mean = ifelse(cl <= 2L, 3, 0), sd = 0.5)
  b <- rnorm(n, mean = ifelse(cl <= 2L, 0.2, 0), sd = 1)
  a[1] <- 8 # extreme selected cell on marker a
  mat <- rbind(b = b, a = a)
  colnames(mat) <- paste0("c", seq_len(n))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(exprs = mat),
    colData = S4Vectors::DataFrame(
      sample_id = rep(c("s1", "s2"), length.out = n),
      cluster_id = factor(cl)
    )
  )
  codes <- t(vapply(seq_len(4), function(k) {
    rowMeans(mat[, cl == k, drop = FALSE])
  }, numeric(2)))
  rownames(codes) <- as.character(seq_len(4))
  S4Vectors::metadata(sce)$SOM_codes <- codes
  S4Vectors::metadata(sce)$map <- list(colsUsed = c("b", "a"))
  sce
}

#' Suppress the expected "plotly event not registered" warning.
#'
#' shiny::testServer() has no live browser client, so plotly's
#' event_register() calls always warn that nothing subscribed to
#' plotly_selected/plotly_relayout events. This is expected noise in
#' headless tests and is suppressed here rather than at each call site.
suppress_plotly_event_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("event tied a source ID.*is not registered", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


# add near the top of the test file, or in helper-*.R
quiet_plotly_test <- function(expr) {
  withCallingHandlers(
    suppressMessages(force(expr)),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("is not registered", msg) ||
        grepl("event_register", msg) ||
        grepl("No trace type specified", msg, ignore.case = TRUE) ||
        grepl("no positional attributes specified", msg, ignore.case = TRUE) ||
        grepl("don't have these attributes", msg, ignore.case = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
