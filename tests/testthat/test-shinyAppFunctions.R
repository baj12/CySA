# Tests for R/shinyAppFunctions.R ----

# ══════════════════════════════════════════════════════════════════════════════
# event wrappers
# ══════════════════════════════════════════════════════════════════════════════

test_that("safe_event_data is exported and aliases exist", {
  expect_type(CySA:::safe_event_data, "closure")
  expect_type(CySA:::.safeEventData, "closure")
})

test_that("quiet_ggplotly returns a plotly object and suppresses harmless warnings", {
  skip_if_not_installed("plotly")
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
    ggplot2::geom_point()
  pl <- quiet_ggplotly(p)
  expect_true(inherits(pl, "plotly") || inherits(pl, "htmlwidget"))
})

test_that("quiet_ggplotly passes source argument to ggplotly", {
  skip_if_not_installed("plotly")
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, customdata = gear)) +
    ggplot2::geom_point()
  pl <- quiet_ggplotly(p, source = "mySource")
  expect_s3_class(pl, "plotly")
})

# ══════════════════════════════════════════════════════════════════════════════
# highlight_df
# ══════════════════════════════════════════════════════════════════════════════

test_that("highlight_df requires metaD", {
  expect_error(highlight_df("x", "y", rs = 1L), "metaD")
})

test_that("highlight_df returns expected columns", {
  metaD <- list(SOM_codes = matrix(
    1:20,
    nrow = 5, ncol = 4,
    dimnames = list(NULL, c("x", "y", "a", "b"))
  ))
  df <- highlight_df("x", "y", rs = c(1L, 3L), metaD = metaD)
  expect_named(df, c("x", "y", "id"))
  expect_equal(df$id, c(1L, 3L))
})


test_that("highlight_df returns empty data.frame for NULL rs", {
  metaD <- list(SOM_codes = matrix(
    1:20,
    nrow = 5, ncol = 4,
    dimnames = list(NULL, c("x", "y", "a", "b"))
  ))
  df <- highlight_df("x", "y", rs = NULL, metaD = metaD)
  expect_equal(nrow(df), 0)
  expect_named(df, c("x", "y", "id"))
})


test_that("highlight_df returns empty data.frame for empty rs", {
  metaD <- list(SOM_codes = matrix(
    1:20,
    nrow = 5, ncol = 4,
    dimnames = list(NULL, c("x", "y", "a", "b"))
  ))
  df <- highlight_df("x", "y", rs = integer(0), metaD = metaD)
  expect_equal(nrow(df), 0)
})


test_that("highlight_df filters out NA and out-of-bounds indices", {
  metaD <- list(SOM_codes = matrix(
    1:20,
    nrow = 5, ncol = 4,
    dimnames = list(NULL, c("x", "y", "a", "b"))
  ))
  df <- highlight_df("x", "y", rs = c(NA, -1, 0, 3, 100), metaD = metaD)
  expect_equal(nrow(df), 1)
  expect_equal(df$id, c(3L))
})


test_that("highlight_df uses somCodesName argument", {
  # Matrix is filled column-wise: [10, 12] in col 1, [11, 13] in col 2
  metaD <- list(
    SOM_codes = matrix(1:4, nrow = 2, dimnames = list(NULL, c("x", "y"))),
    custom_codes = matrix(10:13, nrow = 2, dimnames = list(NULL, c("x", "y")))
  )
  df <- highlight_df("x", "y", rs = c(1L), somCodesName = "custom_codes", metaD = metaD)
  # Row 1, col x = 10; Row 1, col y = 12 (column-wise filling)
  expect_equal(df$x, 10)
  expect_equal(df$y, 12)
})

# ══════════════════════════════════════════════════════════════════════════════
# countBarPlotFunc
# ══════════════════════════════════════════════════════════════════════════════

test_that("countBarPlotFunc returns ggplot for counts-only", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  cpt <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id
  )
  p <- countBarPlotFunc(
    rs = c(1L, 2L),
    clusterPatientTable = cpt,
    cst = "none",
    sce = sce,
    outputList = list(),
    groupsInput = list()
  )
  expect_s3_class(p, "ggplot")
})

test_that("countBarPlotFunc compares to experiment_info numeric column", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  S4Vectors::metadata(sce)$experiment_info$some_numeric <- seq_len(
    nrow(S4Vectors::metadata(sce)$experiment_info)
  )
  cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
  p <- countBarPlotFunc(
    rs = c(1L, 2L), clusterPatientTable = cpt, cst = "some_numeric",
    sce = sce, outputList = list(), groupsInput = list()
  )
  expect_s3_class(p, "ggplot")
})

test_that("countBarPlotFunc compares to named outputList group", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
  p <- countBarPlotFunc(
    rs = c(1L, 2L), clusterPatientTable = cpt, cst = "Ref",
    sce = sce, outputList = list(Ref = c(3L, 4L)), groupsInput = list()
  )
  expect_s3_class(p, "ggplot")
})

test_that("countBarPlotFunc ignores stale node ids in the comparison group", {
  # Regression: a saved group with an id beyond the current SOM's node
  # range crashed with "subscript out of bounds" instead of being ignored.
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
  p <- expect_no_error(countBarPlotFunc(
    rs = c(1L, 2L), clusterPatientTable = cpt, cst = "Ref",
    sce = sce, outputList = list(Ref = c(3L, 4L, 99999L)),
    groupsInput = list()
  ))
  expect_s3_class(p, "ggplot")
})

test_that("countBarPlotFunc errors for unknown comparison column", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
  expect_error(
    countBarPlotFunc(
      rs = c(1L), clusterPatientTable = cpt, cst = "not_found",
      sce = sce, outputList = list(), groupsInput = list()
    ),
    "not found"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# PercentBarPlotFunc
# ══════════════════════════════════════════════════════════════════════════════

test_that("PercentBarPlotFunc returns ggplot", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  S4Vectors::metadata(sce)$experiment_info$some_numeric <- seq_len(
    nrow(S4Vectors::metadata(sce)$experiment_info)
  )
  cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
  p <- PercentBarPlotFunc(
    sce = sce, relativeToCol = "none", clusterPatientTable = cpt,
    rs = c(1L, 2L), outputList = list(), group = NULL,
    groupsInput = list(group1 = character(), group2 = character())
  )
  expect_s3_class(p, "ggplot")
})

test_that("PercentBarPlotFunc respects group colouring", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  S4Vectors::metadata(sce)$experiment_info$group <- rep(
    c("A", "B"),
    length.out = nrow(S4Vectors::metadata(sce)$experiment_info)
  )
  S4Vectors::metadata(sce)$experiment_info$some_numeric <- seq_len(
    nrow(S4Vectors::metadata(sce)$experiment_info)
  )
  cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)

  expect_no_warning(
    p <- PercentBarPlotFunc(
      sce = sce, relativeToCol = "none", clusterPatientTable = cpt,
      rs = c(1L, 2L), outputList = list(), group = NULL,
      groupsInput = list(group1 = "sample1", group2 = "sample2")
    )
  )
  expect_s3_class(p, "ggplot")
})

# ══════════════════════════════════════════════════════════════════════════════
# ggsomPlot / somPlot
# ══════════════════════════════════════════════════════════════════════════════

test_that("ggsomPlot returns ggplot with dimSelection", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  metaD <- S4Vectors::metadata(sce)
  pp1 <- plotSOMScatter(sce, chs = c("marker1", "marker2"))
  dimSelection <- list(list(
    dims = c("marker1", "marker2"),
    xlim = c(0, 1), ylim = c(0, 1),
    xzoom = NULL, yzoom = NULL
  ))
  p <- ggsomPlot(
    pp1,
    plotIdx = 1L, rs = c(1L, 2L),
    dimSelection = dimSelection, sce = sce, metaD = metaD
  )
  expect_s3_class(p, "ggplot")
})

test_that("ggsomPlot applies explicit xlim/ylim", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  metaD <- S4Vectors::metadata(sce)
  pp1 <- plotSOMScatter(sce, chs = c("marker1", "marker2"))
  dimSelection <- list(list(
    dims = c("marker1", "marker2"),
    xlim = c(0, 1), ylim = c(0, 1),
    xzoom = NULL, yzoom = NULL
  ))
  p <- ggsomPlot(
    pp1,
    plotIdx = 1L, rs = c(1L),
    dimSelection = dimSelection, sce = sce, metaD = metaD,
    xlim = c(-1, 2), ylim = c(-1, 2)
  )
  expect_s3_class(p, "ggplot")
})

test_that("somPlot returns NULL when pp1 is NULL", {
  expect_null(somPlot(
    pp1 = NULL, plotIdx = 1L, rs = 1L, colorbyGroups = "",
    showGroups = FALSE, sce = CySA_example_sce(n_cells = 100, n_nodes = 4)
  ))
})

test_that("somPlot returns ggplotly-compatible object", {
  skip_if_not_installed("plotly")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  pp1 <- plotSOMScatter(sce, chs = c("marker1", "marker2"))
  dimSelection <- list(list(
    dims = c("marker1", "marker2"),
    xlim = c(0, 1), ylim = c(0, 1),
    xzoom = NULL, yzoom = NULL
  ))
  p <- somPlot(
    pp1,
    plotIdx = 1L, rs = c(1L, 2L),
    colorbyGroups = NULL, showGroups = FALSE,
    dimSelection = dimSelection, sce = sce,
    somCodesName = "SOM_codes"
  )
  expect_true(inherits(p, "plotly") || inherits(p, "ggplot"))
})

test_that(".somScatterLegendOn fires for known-type colours only", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  # Numeric statistics and node dispersion carry a colour bar, not a legend.
  expect_false(.somScatterLegendOn(sce, "n"))
  expect_false(.somScatterLegendOn(sce, "node_dispersion"))
  labeled <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", rep(1:3, length.out = 10)))
  labeled <- prepClusterSelectorData(labeled, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")$sce
  expect_true(.somScatterLegendOn(labeled, "cell_type_majority"))
  expect_false(.somScatterLegendOn(labeled, "purity"))
})

test_that("somPlot shows the legend for known-type colours only", {
  skip_if_not_installed("plotly")
  labeled <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", rep(1:3, length.out = 10)))
  labeled <- prepClusterSelectorData(labeled, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")$sce
  dimSelection <- list(list(
    dims = c("marker1", "marker2"),
    xlim = c(0, 1), ylim = c(0, 1),
    xzoom = NULL, yzoom = NULL
  ))
  pp1 <- plotSOMScatter(labeled, chs = c("marker1", "marker2"),
    color_by = "cell_type_majority")
  p <- somPlot(pp1, plotIdx = 1L, rs = c(1L, 2L),
    colorbyGroups = NULL, showGroups = FALSE, dimSelection = dimSelection,
    sce = labeled, somCodesName = "SOM_codes",
    colorVar = "cell_type_majority")
  expect_true(isTRUE(p$x$layout$showlegend))
  # Numeric colourings keep the legend hidden.
  pp2 <- plotSOMScatter(labeled, chs = c("marker1", "marker2"))
  p2 <- somPlot(pp2, plotIdx = 1L, rs = c(1L, 2L),
    colorbyGroups = NULL, showGroups = FALSE, dimSelection = dimSelection,
    sce = labeled, somCodesName = "SOM_codes", colorVar = "n")
  expect_false(isTRUE(p2$x$layout$showlegend))
})

# ══════════════════════════════════════════════════════════════════════════════
# tsneFunc
# ══════════════════════════════════════════════════════════════════════════════

test_that("tsneFunc returns Rtsne object", {
  skip_if_not_installed("Rtsne")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  res <- tsneFunc(
    dimRedSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:4],
    perplexity = 2L, sce = sce
  )
  expect_s3_class(res, "Rtsne")
})

test_that("tsneFunc clamps perplexity to max allowed", {
  skip_if_not_installed("Rtsne")
  sce <- CySA_example_sce(n_cells = 100, n_nodes = 4)
  res <- tsneFunc(
    dimRedSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:4],
    perplexity = 999L, sce = sce
  )
  expect_s3_class(res, "Rtsne")
})

# ══════════════════════════════════════════════════════════════════════════════
# compute_relative_counts
# ══════════════════════════════════════════════════════════════════════════════

test_that("compute_relative_counts handles 'none' mode", {
  cpt <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("1", "2"))
  )
  class(cpt) <- "table"
  res <- compute_relative_counts(cpt,
    rs = "1", relativeToCol = "none",
    expInfo = data.frame(), outputList = list()
  )
  expect_type(res, "double")
  expect_length(res, 2L)
})

test_that("compute_relative_counts divides by experiment_info column", {
  cpt <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("1", "2"))
  )
  class(cpt) <- "table"
  expInfo <- data.frame(total = c(100, 200), row.names = c("s1", "s2"))
  res <- compute_relative_counts(cpt,
    rs = "1", relativeToCol = "total",
    expInfo = expInfo, outputList = list()
  )
  expect_equal(unname(res), c(10, 10))
})

test_that("compute_relative_counts divides by outputList group", {
  cpt <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("1", "2"))
  )
  class(cpt) <- "table"
  res <- compute_relative_counts(cpt,
    rs = "1", relativeToCol = "Ref",
    expInfo = data.frame(),
    outputList = list(Ref = c("2"))
  )
  expect_equal(unname(res), c(10 / 30, 20 / 40) * 100)
})

test_that("compute_relative_counts errors for unknown relativeToCol", {
  cpt <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"), c("1", "2")))
  class(cpt) <- "table"
  expect_error(
    compute_relative_counts(cpt,
      rs = "1", relativeToCol = "missing",
      expInfo = data.frame(), outputList = list()
    ),
    "not found"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# buildProjectionDf / drawProjection
# ══════════════════════════════════════════════════════════════════════════════

test_that("buildProjectionDf requires som_codes and SOM_stats", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  S4Vectors::metadata(sce)$SOM_stats <- NULL
  dimSelection <- list(list(dims = c("marker1", "marker2")))
  # shiny::req() stops silently; verify the function does not return a data frame.
  expect_warning(
    expect_error(
      buildProjectionDf(ggplot2::ggplot(), plotIdx = 1L, dimSelection = dimSelection, sce = sce),
      regexp = NA
    ), "not found in metadata"
  )
})

test_that("buildProjectionDf returns expected column order", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  dimSelection <- list(list(dims = c("marker1", "marker2")))
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, dimSelection, sce)
  expect_equal(names(df)[1:2], c("marker1", "marker2"))
  expect_true("id" %in% names(df))
})

test_that("buildProjectionDf joins node_labels and preserves NA majority", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", rep(1:2, 5)))
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  lab <- prepped$sce$cell_type
  lab[prepped$sce$cluster_id %in% 9:10] <- NA
  prepped$sce$cell_type <- lab
  prepped$sce <- .applyLabelStats(prepped$sce, "cell_type")
  dimSelection <- list(list(dims = c("marker1", "marker2")))
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, dimSelection, prepped$sce)
  expect_true(all(c("purity", "cell_type_majority") %in% colnames(df)))
  expect_true(any(is.na(df$cell_type_majority[df$id %in% 9:10])))
  # Unknown purity must not be zeroed into a fake 0% pure node.
  expect_true(all(is.na(df$purity[df$id %in% 9:10])))
  expect_true(all(is.finite(df$purity[df$id %in% 1:8])))
})

test_that("buildProjectionDf numeric NA-to-0 behavior is unchanged", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  dimSelection <- list(list(dims = c("marker1", "marker2")))
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, dimSelection, sce)
  expect_false(anyNA(df))
})

test_that("drawProjection returns ggplot without colorbyGroups", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, list(list(dims = c("marker1", "marker2"))), sce)

  p <- drawProjection(df, rs = c(1L, 2L), colorbyGroups = NULL, sce = sce)
  expect_s3_class(p, "ggplot")
})

test_that("drawProjection colours by groups when colorbyGroups supplied", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, list(list(dims = c("marker1", "marker2"))), sce)

  p <- drawProjection(
    df,
    rs = 1L, colorbyGroups = c("A"),
    sce = sce, outputList = list(A = 1:2, B = 3:4)
  )

  expect_s3_class(p, "ggplot")
})

test_that("drawProjection colours by agreement when supplied", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  codes <- S4Vectors::metadata(sce)$SOM_codes
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, list(list(dims = c("marker1", "marker2"))), sce)
  xy <- as.matrix(df[, c("marker1", "marker2")])
  agreement <- .knnAgreement(codes, xy, k = 2)

  p <- drawProjection(df, rs = 1L, colorbyGroups = NULL, sce = sce,
    agreement = agreement)

  expect_s3_class(p, "ggplot")
  # The colour legend title must name the score; labels materialise at
  # build time.
  b <- ggplot2::ggplot_build(p)
  expect_identical(b$plot$labels$colour, "agreement")
})

test_that("drawProjection group colours win over agreement", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  codes <- S4Vectors::metadata(sce)$SOM_codes
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, list(list(dims = c("marker1", "marker2"))), sce)
  xy <- as.matrix(df[, c("marker1", "marker2")])
  agreement <- .knnAgreement(codes, xy, k = 2)

  p <- drawProjection(df, rs = 1L, colorbyGroups = c("A"),
    sce = sce, outputList = list(A = 1:2, B = 3:4), agreement = agreement)

  # Explicit group colouring must not be overridden by agreement mode:
  # the colour legend stays the group variable, not the score.
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_false(grepl("agreement", b$plot$labels$colour %||% "",
    ignore.case = TRUE))
})

test_that("drawProjection agreement tooltip adds the score line", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  codes <- S4Vectors::metadata(sce)$SOM_codes
  df <- buildProjectionDf(ggplot2::ggplot(), 1L, list(list(dims = c("marker1", "marker2"))), sce)
  xy <- as.matrix(df[, c("marker1", "marker2")])
  agreement <- .knnAgreement(codes, xy, k = 2)

  p <- drawProjection(df, rs = 1L, colorbyGroups = NULL, sce = sce,
    agreement = agreement)
  b <- ggplot2::ggplot_build(p)$data[[1]]
  # The tooltip gains an "agreement:" line with the score in [0, 1]
  # (ggplot2 consumes the raw column into the colour mapping, so the
  # tooltip text is the observable).
  expect_true(all(grepl("agreement:", b$text, fixed = TRUE)))
  scores <- as.numeric(sub(".*agreement: ([0-9.]+).*", "\\1", b$text))
  expect_true(all(scores >= 0 & scores <= 1))
})

# ══════════════════════════════════════════════════════════════════════════════
# .knnAgreement: neighbourhood agreement between a 2D view and full codes
# ══════════════════════════════════════════════════════════════════════════════

test_that(".knnAgreement scores neighbourhood overlap per node", {
  # Four nodes: codes place 1 close to 2 and 3 close to 4 in full space.
  codes <- rbind(c(0.0, 0.0), c(0.1, 0.0), c(5.0, 5.0), c(5.1, 5.0))
  # The 2D view keeps the pair structure perfectly: same neighbours.
  xy_good <- rbind(c(0.0, 0.0), c(0.1, 0.0), c(5.0, 5.0), c(5.1, 5.0))
  scores_good <- .knnAgreement(codes, xy_good, k = 1)
  expect_equal(unname(scores_good), c(1, 1, 1, 1))
  # The 2D view destroys the pair structure: every node's nearest
  # neighbour is wrong, so all scores are 0.
  xy_bad <- rbind(c(0.0, 0.0), c(5.0, 5.0), c(5.1, 5.0), c(0.1, 0.0))
  scores_bad <- .knnAgreement(codes, xy_bad, k = 1)
  expect_equal(unname(scores_bad), c(0, 0, 0, 0))
})

test_that(".knnAgreement with k >= n returns agreement 1 everywhere", {
  codes <- rbind(c(0, 0), c(0.1, 0), c(5, 5))
  xy <- rbind(c(3, 3), c(3, 3.1), c(3, 3.2))
  scores <- .knnAgreement(codes, xy, k = 3)
  expect_equal(unname(scores), c(1, 1, 1))
})

test_that(".knnAgreement handles ties deterministically", {
  # Three nodes mutually equidistant: with k = 1, each node's kNN set in
  # both spaces must contain the same single neighbour regardless of how
  # ties break, so agreement is 1 everywhere.
  codes <- rbind(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2))
  xy <- rbind(c(10, 10), c(11, 10), c(10.5, 10 + sqrt(3) / 2))
  scores <- .knnAgreement(codes, xy, k = 1)
  expect_equal(unname(scores), c(1, 1, 1))
})

test_that(".knnAgreement drops NA coordinates instead of failing", {
  codes <- rbind(c(0.0, 0.0), c(0.1, 0.0), c(5.0, 5.0), c(5.1, 5.0))
  xy <- rbind(c(0.0, 0.0), c(0.1, 0.0), c(5.0, 5.0), c(5.1, 5.0))
  xy[2, 1] <- NA
  scores <- .knnAgreement(codes, xy, k = 1)
  expect_true(is.na(scores[2]))
  expect_equal(unname(scores[c(1, 3, 4)]), c(1, 1, 1))
})

test_that(".knnAgreement is 1 for a single node", {
  # No neighbours exist in either space: the degenerate score is 1
  # (nothing to distort) rather than NaN.
  scores <- .knnAgreement(rbind(c(0, 0)), rbind(c(1, 1)), k = 5)
  expect_equal(unname(scores), 1)
})

# ══════════════════════════════════════════════════════════════════════════════
# violin plots
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotViolinFunc returns ggplot with enough data", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 4, n_markers = 6)
  outputList <- list(A = c(1L, 2L), B = c(3L, 4L))
  p <- plotViolinFunc(
    sce = sce, upsetSelection = names(outputList),
    outputList = outputList,
    violinSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:3]
  )
  expect_s3_class(p, "ggplot")
})

test_that("plotViolinFunc returns placeholder when upsetSelection too short", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  outputList <- list(A = c(1L, 2L), B = c(3L, 4L))
  expect_message(
    p <- plotViolinFunc(
      sce = sce, upsetSelection = "A", outputList = outputList,
      violinSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:2]
    ),
    "please check that all upsetSelection are in outputList:"
  )
  expect_true(is.null(p) || inherits(p, "ggplot"))
})

test_that("plotViolinFunc returns empty plot when no groups have >=2 nodes", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  outputList <- list(A = 1L)
  p <- plotViolinFunc(
    sce = sce, upsetSelection = "A", outputList = outputList,
    violinSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:2]
  )
  expect_s3_class(p, "ggplot")
})

test_that("plotViolinFunc ignores a stale node id in a saved group", {
  # Regression: somCodes[rows, ...] with a row id beyond nrow(somCodes)
  # threw "subscript out of bounds" instead of dropping the stale id.
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 4, n_markers = 6)
  outputList <- list(A = c(1L, 2L), B = c(3L, 4L, 99999L))
  p <- expect_no_error(plotViolinFunc(
    sce = sce, upsetSelection = names(outputList),
    outputList = outputList,
    violinSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:3]
  ))
  expect_s3_class(p, "ggplot")
})

test_that("plotViolin2Func returns ggplot", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 4, n_markers = 6)
  outputList <- list(A = c(1L, 2L), B = c(3L, 4L))
  p <- plotViolin2Func(
    sce = sce,
    violinSelection = colnames(S4Vectors::metadata(sce)$SOM_codes)[1:3],
    upsetSelection = names(outputList),
    outputList = outputList
  )
  expect_s3_class(p, "ggplot")
})

# ══════════════════════════════════════════════════════════════════════════════
# upsetPlotFunc
# ══════════════════════════════════════════════════════════════════════════════

test_that("upsetPlotFunc returns NULL for fewer than two groups", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if(
    R.version$status == "Under development (unstable)",
    paste(
      "Known ComplexHeatmap/circlize 'node stack overflow' on R-devel",
      "(reproduces with ComplexHeatmap 2.29.0 + circlize 0.4.18 on both",
      "linux-devel and windows-devel CI runners; passes on release/oldrel).",
      "Root cause appears tied to degenerate/tied input values feeding",
      "ComplexHeatmap's internal dendrogram/ordering logic -- see",
      "jokergoo/ComplexHeatmap issues #199, #539, #541 for the historical",
      "pattern. Not reproducible locally on R-release."
    )
  )
  expect_null(upsetPlotFunc("A", list(A = 1L), CySA_example_sce(n_cells = 100, n_nodes = 4)))
})

test_that("upsetPlotFunc returns ComplexHeatmap object for valid groups", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if(
    R.version$status == "Under development (unstable)",
    paste(
      "Known ComplexHeatmap/circlize 'node stack overflow' on R-devel",
      "(reproduces with ComplexHeatmap 2.29.0 + circlize 0.4.18 on both",
      "linux-devel and windows-devel CI runners; passes on release/oldrel).",
      "Root cause appears tied to degenerate/tied input values feeding",
      "ComplexHeatmap's internal dendrogram/ordering logic -- see",
      "jokergoo/ComplexHeatmap issues #199, #539, #541 for the historical",
      "pattern. Not reproducible locally on R-release."
    )
  )

  sce <- CySA_example_sce(n_cells = 300, n_nodes = 6)
  S4Vectors::metadata(sce)$SOM_stats$n <- c(37L, 61L, 45L, 58L, 29L, 70L)
  outputList <- list(A = c(1L, 2L, 3L), B = c(2L, 3L, 5L))

  p <- upsetPlotFunc(names(outputList), outputList, sce)
  expect_s4_class(p, "Heatmap")
})
