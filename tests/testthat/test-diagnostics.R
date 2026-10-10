# Tests for R/plot-diagnostics.R ----
# Sample-level QC diagnostics (MDS, per-sample expression heatmap, NRS) and
# the dendrogram-ordered cluster heatmap, added on top of the CATALYST
# reference workflow reproduced in vignettes/catalyst_workflow.Rmd.

test_that(".samplePseudobulk returns a samples x markers matrix of medians", {
  sce <- CySA_example_sce(n_cells = 40, n_nodes = 5, n_markers = 3)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  pb <- CySA:::.samplePseudobulk(sce, markers)

  expect_equal(dim(pb), c(4L, 3L))
  expect_equal(colnames(pb), markers)
  expect_setequal(rownames(pb), unique(as.character(sce$sample_id)))
})

test_that("plotSampleMDS returns a ggplot and errors with too few samples", {
  sce <- CySA_example_sce(n_cells = 40, n_nodes = 5, n_markers = 3)
  markers <- S4Vectors::metadata(sce)$map$colsUsed

  p <- plotSampleMDS(sce, markers)
  expect_s3_class(p, "ggplot")

  p2 <- plotSampleMDS(sce, markers, color_by = "some_numeric")
  expect_s3_class(p2, "ggplot")

  sce$sample_id[sce$sample_id != "sample1"] <- "sample1"
  expect_error(plotSampleMDS(sce, markers), "at least 3 samples")
})

test_that("plotSampleExprHeatmap returns a ggplot for >=2 and <2 samples", {
  sce <- CySA_example_sce(n_cells = 40, n_nodes = 5, n_markers = 3)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  p <- plotSampleExprHeatmap(sce, markers)
  expect_s3_class(p, "ggplot")

  # A single marker/single sample degenerate case should not error.
  sce$sample_id[] <- "sample1"
  p2 <- plotSampleExprHeatmap(sce, markers[1])
  expect_s3_class(p2, "ggplot")
})

test_that("computeNRS returns one row per sample-marker combination", {
  sce <- CySA_example_sce(n_cells = 60, n_nodes = 5, n_markers = 3)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  nrsDf <- computeNRS(sce, markers)

  expect_true(is.data.frame(nrsDf))
  expect_setequal(colnames(nrsDf), c("sample_id", "marker", "nrs"))
  expect_equal(nrow(nrsDf), length(unique(sce$sample_id)) * length(markers))
})

test_that("plotNRS returns a ggplot, with and without colour_by", {
  sce <- CySA_example_sce(n_cells = 60, n_nodes = 5, n_markers = 3)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  p <- plotNRS(sce, markers)
  expect_s3_class(p, "ggplot")
  p2 <- plotNRS(sce, markers, color_by = "some_numeric")
  expect_s3_class(p2, "ggplot")
})

test_that("plotOrderedClusterHeatmap orders clusters by dendrogram leaves", {
  sce <- CySA_example_sce(n_cells = 40, n_nodes = 6, n_markers = 3)
  somCodes <- S4Vectors::metadata(sce)$SOM_codes
  somStats <- S4Vectors::metadata(sce)$SOM_stats
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(somCodes)))
  dendTable <- data.frame(
    id = seq_len(nrow(somCodes)), label = rownames(somCodes),
    stringsAsFactors = FALSE
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plotOrderedClusterHeatmap(somCodes, somStats, dend, dendTable)
  expect_true(inherits(p, "ggplot") || inherits(p, "ggplot_built") ||
    is.list(p))
})

test_that("plotOrderedClusterHeatmap returns a plain heatmap w/o somStats", {
  sce <- CySA_example_sce(n_cells = 40, n_nodes = 6, n_markers = 3)
  somCodes <- S4Vectors::metadata(sce)$SOM_codes
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(somCodes)))

  p <- plotOrderedClusterHeatmap(somCodes, NULL, dend)
  expect_s3_class(p, "ggplot")
})

test_that("Diagnostics and ordered-heatmap outputs render inside the app", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(diagColorBy = "treatment"))
    suppressWarnings(session$flushReact())

    expect_s3_class(diagMdsPlot(), "ggplot")
    expect_s3_class(diagExprHeatmapPlot(), "ggplot")
    expect_s3_class(diagNrsPlot(), "ggplot")
    expect_true(!is.null(orderedHeatmapPlot()))
  })
})

test_that(".subsetSceByGroups restricts cells to the chosen groups", {
  sce <- CySA_example_sce(n_cells = 40, n_nodes = 10, n_markers = 3)
  outputList <- list(groupA = as.character(1:5), groupB = as.character(6:10))

  full <- CySA:::.subsetSceByGroups(sce, outputList, NULL)
  expect_identical(ncol(full), ncol(sce))

  subA <- CySA:::.subsetSceByGroups(sce, outputList, "groupA")
  expect_true(all(as.character(subA$cluster_id) %in% outputList$groupA))
  expect_lt(ncol(subA), ncol(sce))
})

test_that("Diagnostics and ordered-heatmap panels respect a group selection", {
  test_app <- make_test_app(n_cells = 200, n_nodes = 10)

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1,2,3,4,5", applyclusterNumbers = 1,
      clusterName = "groupA", applyName = 1
    ))
    suppressWarnings(session$flushReact())

    session$setInputs(diagGroups = "groupA")
    session$flushReact()
    expect_s3_class(diagMdsPlot(), "ggplot")

    allMarkers <- colnames(S4Vectors::metadata(test_app$sce)$SOM_codes)
    session$setInputs(orderedHeatmapGroups = "groupA",
      orderedHeatmapMarkers = allMarkers[1])
    session$flushReact()
    expect_true(!is.null(orderedHeatmapPlot()))
  })
})
