# Additional coverage tests for plotting functions ----

# ══════════════════════════════════════════════════════════════════════════════
# plotCytoScatter
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotCytoScatter respects zeros = TRUE", {
  skip_if_not_installed("KernSmooth")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  p <- plotCytoScatter(sce, chs = c("marker1", "marker2"), zeros = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plotCytoScatter colours by numeric colData column", {
  skip_if_not_installed("KernSmooth")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  sce$num_col <- rnorm(ncol(sce))
  p <- plotCytoScatter(sce, chs = c("marker1", "marker2"), color_by = "num_col")
  expect_s3_class(p, "ggplot")
})

test_that("plotCytoScatter colours by cluster column via CATALYST", {
  skip_if_not_installed("KernSmooth")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  p <- plotCytoScatter(sce, chs = c("marker1", "marker2"), color_by = "cluster_id")
  expect_s3_class(p, "ggplot")
})

test_that("plotCytoScatter facets when >2 channels", {
  skip_if_not_installed("KernSmooth")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  p <- plotCytoScatter(sce, chs = c("marker1", "marker2", "marker3", "marker4"))
  expect_s3_class(p, "ggplot")
})

test_that("plotCytoScatter uses label = 'target'", {
  skip_if_not_installed("KernSmooth")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  p <- plotCytoScatter(sce, chs = c("marker1", "marker2"), label = "target")
  expect_s3_class(p, "ggplot")
})

# ══════════════════════════════════════════════════════════════════════════════
# plotSOMScatter
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotSOMScatter warns when statsSlot is missing", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  S4Vectors::metadata(sce)$SOM_stats <- NULL
  expect_warning(
    plotSOMScatter(sce, chs = c("marker1", "marker2"), statsSlot = "SOM_stats"),
    "SOM_stats"
  )
})

test_that("plotSOMScatter computes missing SOM code columns", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4, n_markers = 8)
  keep <- colnames(S4Vectors::metadata(sce)$SOM_codes)[1:2]
  S4Vectors::metadata(sce)$SOM_codes <- S4Vectors::metadata(sce)$SOM_codes[, keep, drop = FALSE]
  expect_warning(
    expect_warning(
      p <- plotSOMScatter(sce, chs = c("marker3", "marker4")),
      "computing SOM stats"
    ),
    "computing SOM stats"
  )
  expect_s3_class(p, "ggplot")
})

test_that("plotSOMScatter facets when >2 channels", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  p <- plotSOMScatter(
    sce,
    chs = c("marker1", "marker2", "marker3", "marker4")
  )
  expect_s3_class(p, "ggplot")
})

test_that("plotSOMScatter returns NULL when all channels unknown", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  expect_warning(
    p <- plotSOMScatter(sce, chs = c("no", "way")),
    "Unknown channels/markers"
  )
  expect_null(p)
})

# ══════════════════════════════════════════════════════════════════════════════
# plotSOMScatter with known-label node properties
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotSOMScatter colors by purity and by known cell type", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", rep(1:3, length.out = 10)))
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  get_scale <- function(p, aesthetic) {
    Filter(function(s) aesthetic %in% s$aesthetics, p$scales$scales)[[1]]
  }
  p1 <- plotSOMScatter(prepped$sce, chs = c("marker1", "marker2"),
    color_by = "purity")
  # ggplot2 4.x names colour scales ScaleContinuous/ScaleDiscrete (no
  # *Colour subclass); the aesthetic filter pins them to colour.
  expect_s3_class(get_scale(p1, "colour"), "ScaleContinuous")
  p2 <- plotSOMScatter(prepped$sce, chs = c("marker1", "marker2"),
    color_by = "cell_type_majority")
  expect_s3_class(get_scale(p2, "colour"), "ScaleDiscrete")
})

test_that("plotSOMScatter ignores unknown color_by without error", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  expect_s3_class(plotSOMScatter(sce, chs = c("marker1", "marker2"),
    color_by = "purity"), "ggplot")
})

test_that("NA-purity nodes stay plotted when purity is the size", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", rep(1:3, length.out = 10)))
  # Node 3 holds only NA-labelled cells: its purity is NA.
  lab <- sce$cell_type
  lab[sce$cluster_id == "3"] <- NA
  sce$cell_type <- lab
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  nl <- S4Vectors::metadata(prepped$sce)$node_labels
  expect_true(is.na(nl$purity[nl$id == 3L]))
  p <- plotSOMScatter(prepped$sce, chs = c("marker1", "marker2"),
    pointSize = "purity", color_by = "purity", zeros = TRUE)
  built <- ggplot2::ggplot_build(p)
  stats <- S4Vectors::metadata(prepped$sce)$SOM_stats
  # Every node is drawn, with a finite mapped size for the NA-purity node.
  expect_identical(nrow(built$data[[1]]), nrow(stats))
  expect_true(all(is.finite(built$data[[1]]$size)))
})

# ══════════════════════════════════════════════════════════════════════════════
# .somStatChoices (sidebar color/size choices)
# ══════════════════════════════════════════════════════════════════════════════

test_that("somStatChoices exposes label columns only when labels exist", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  ch <- .somStatChoices(sce)
  expect_identical(ch$colorStats, c("n", "mean", "median", "rdQu", "max"))
  expect_identical(ch$sizeStats, c("n", "mean", "median", "rdQu", "max"))
  sce2 <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", 1:10))
  prepped <- prepClusterSelectorData(sce2, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  ch2 <- .somStatChoices(prepped$sce)
  expect_setequal(ch2$colorStats,
    c("n", "mean", "median", "rdQu", "max", "purity", "cell_type_majority"))
  expect_true("purity" %in% ch2$sizeStats)
  expect_false("cell_type_majority" %in% ch2$sizeStats)
})

# ══════════════════════════════════════════════════════════════════════════════
# .buildNodePropsBox / .nodePropertyChoices (node-properties box)
# ══════════════════════════════════════════════════════════════════════════════

test_that("node properties box renders and prunes label choices", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  box1 <- .buildNodePropsBox(.nodePropertyChoices(S4Vectors::metadata(sce)),
    hasLabels = FALSE)
  expect_s3_class(box1, "shiny.tag")
  txt <- paste(capture.output(print(box1)), collapse = "")
  expect_true(grepl("nodePropsX", txt))
  expect_true(grepl("nodePropsY", txt))
  expect_true(grepl("nodePropsPlot", txt))
  # no purity among the axis choices (tooltips may still mention the word)
  expect_false(grepl("<option value=\"purity\"", txt))
  sce2 <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", 1:10))
  prepped <- prepClusterSelectorData(sce2, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  box2 <- .buildNodePropsBox(
    .nodePropertyChoices(S4Vectors::metadata(prepped$sce)),
    hasLabels = TRUE)
  txt2 <- paste(capture.output(print(box2)), collapse = "")
  expect_true(grepl("<option value=\"purity\"", txt2))
})

# ══════════════════════════════════════════════════════════════════════════════
# validate-inputs
# ══════════════════════════════════════════════════════════════════════════════

test_that(".validateClusterSelectorInputs checks fsom type", {
  skip_if_not_installed("FlowSOM")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  # Create proper somRasterData with all marker columns
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(S4Vectors::metadata(sce)$SOM_codes)),
    y = rep(seq_len(2), each = nrow(S4Vectors::metadata(sce)$SOM_codes) / 2)
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(S4Vectors::metadata(sce)$SOM_codes)) / nrow(S4Vectors::metadata(sce)$SOM_codes)
  }
  expect_error(
    CySA:::.validateClusterSelectorInputs(
      sce, sce, list(), list(markers[1:2]),
      stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
      data.frame(id = 1, label = "node1"),
      table(sample_id = sce$sample_id, cluster_id = sce$cluster_id),
      somRasterData,
      "SOM_codes", 1L,
      fsom = "not_a_fsom"
    ),
    "FlowSOM"
  )
})

test_that(".validateClusterSelectorInputs checks fsom MST layout", {
  skip_if_not_installed("FlowSOM")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  # Create a FlowSOM object with missing MST$l layout
  bad_fsom <- list(MST = list()) # Missing 'l' component
  class(bad_fsom) <- "FlowSOM"
  # Create proper somRasterData with all marker columns
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(S4Vectors::metadata(sce)$SOM_codes)),
    y = rep(seq_len(2), each = nrow(S4Vectors::metadata(sce)$SOM_codes) / 2)
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(S4Vectors::metadata(sce)$SOM_codes)) / nrow(S4Vectors::metadata(sce)$SOM_codes)
  }
  expect_error(
    CySA:::.validateClusterSelectorInputs(
      sce,
      sce_subsampled = sce,
      outputList = list(),
      dList = list(markers[1:2]),
      dend = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
      dendTable = data.frame(id = 1, label = "node1"),
      clusterPatientTable = table(sample_id = sce$sample_id, cluster_id = sce$cluster_id),
      somRasterData = somRasterData,
      somCodesName = "SOM_codes",
      nPlots = 1L, fsom = bad_fsom
    ),
    "MST"
  )
})

test_that(".validateClusterSelectorInputs checks somRasterData columns", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  bad_raster <- data.frame(x = 1, y = 1)
  expect_error(
    CySA:::.validateClusterSelectorInputs(
      sce, sce, list(), list(c("marker1", "marker2")),
      stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
      data.frame(id = 1, label = "node1"),
      table(sample_id = sce$sample_id, cluster_id = sce$cluster_id),
      bad_raster,
      "SOM_codes", 1L,
      fsom = NULL
    ),
    "missing SOM columns"
  )
})
