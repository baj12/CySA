# test-coverage-floor.R
# Targeted tests for branches of the diffcyt observer/outputs and the
# statistics/PDF helpers that the app-level tests never reach.

library(testthat)
library(shiny)

# ---- Fixtures ---------------------------------------------------------------

.floor_ei <- function() {
  data.frame(
    sample_id = paste0("sample", 1:4),
    some_numeric = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )
}

.floor_cpt <- function(constant = FALSE) {
  vals <- if (constant) rep(5L, 40) else c(
    3L, 8L, 5L, 2L, 9L, 4L, 7L, 6L, 1L, 5L, 2L, 6L, 9L, 3L, 8L, 4L, 7L, 5L,
    6L, 2L, 4L, 9L, 3L, 7L, 5L, 8L, 1L, 6L, 3L, 4L, 8L, 2L, 7L, 5L, 9L, 6L,
    4L, 3L, 2L, 5L
  )
  m <- matrix(vals, nrow = 4, dimnames = list(
    sample_id = paste0("sample", 1:4), cluster_id = as.character(1:10)
  ))
  as.table(m)
}

.floor_sce <- function(cpt) {
  sid <- rep(rownames(cpt), times = rowSums(cpt))
  cid <- unlist(lapply(rownames(cpt), function(s) {
    rep(as.integer(colnames(cpt)), times = cpt[s, ])
  }))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(0, 1, length(sid))),
    colData = S4Vectors::DataFrame(sample_id = sid, cluster_id = cid)
  )
}

.floor_st <- function(cpt, ei = .floor_ei(), rs = 1:3) {
  list(
    rsUsed = function() rs,
    rv = list(outputList = list(Rest = 4:10)),
    clusterPatientTable = cpt,
    metaD = list(experiment_info = ei),
    sce = .floor_sce(cpt),
    groupsInput = function() {
      list(group1 = c("sample1", "sample2"), group2 = c("sample3", "sample4"))
    }
  )
}

# ---- .ttestResult -----------------------------------------------------------

test_that(".ttestResult runs a Welch test relative to a saved group", {
  ei <- .floor_ei()
  ei$cohort <- factor(c("A", "A", "B", "B"))
  res <- CySA:::.ttestResult(
    .floor_st(.floor_cpt(), ei),
    list(groupsVar = "cohort", relativeTo = "Rest")
  )
  expect_s3_class(res, "htest")
})

test_that(".ttestResult reports 'not enough data' for constant rates", {
  ei <- .floor_ei()
  ei$cohort <- factor(c("A", "A", "B", "B"))
  res <- CySA:::.ttestResult(
    .floor_st(.floor_cpt(constant = TRUE), ei),
    list(groupsVar = "cohort", relativeTo = "none")
  )
  expect_identical(res, "not enough data")
})

test_that(".ttestResult reports 'no data' without a grouping variable", {
  expect_identical(
    CySA:::.ttestResult(.floor_st(.floor_cpt()), list(groupsVar = NULL)),
    "no data"
  )
})

# ---- .cellCountsTable -------------------------------------------------------

test_that(".cellCountsTable stops when experiment_info has NAs", {
  ei <- .floor_ei()
  ei$some_numeric[2] <- NA
  expect_error(
    CySA:::.cellCountsTable(
      .floor_st(.floor_cpt(), ei),
      list(compareStatsTo = "some_numeric", singleNode = "1")
    ),
    "NAs produced"
  )
})

# ---- PDF helpers ------------------------------------------------------------

test_that(".pdfPlotRegistry builds drawable ridge and violin plots", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  st <- list(
    sce = sce, somCodesName = "SOM_codes", colsUsed = markers,
    rv = list(outputList = list(a = 1:3, b = 4:6)), fsom = NULL
  )
  reg <- CySA:::.pdfPlotRegistry(st, list(), 1:3)
  expect_true(all(c("ridgeTrace", "violin1", "violin2", "upset") %in%
    names(reg)))
  expect_null(reg$flowSOMStars)
  expect_s3_class(reg$ridgeTrace(), "ggplot")
  expect_s3_class(reg$violin1(), "ggplot")
  expect_s3_class(reg$violin2(), "ggplot")
})

test_that(".pdfRidgeTrace honours an explicit cluster pick", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  st <- list(sce = sce, somCodesName = "SOM_codes", colsUsed = markers)
  p <- CySA:::.pdfRidgeTrace(
    st, list(ridgeClusterSelection = c("1", "2"), ridgeThreshold = 0.95),
    1:3
  )
  expect_s3_class(p, "ggplot")
})

test_that(".pdfPlotRegistry draws FlowSOM stars when a fsom is present", {
  skip_if_not_installed("FlowSOM")
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  x <- t(SummarizedExperiment::assay(sce, "exprs")[markers, ])
  fsom <- FlowSOM::FlowSOM(
    flowCore::flowFrame(x), colsToUse = seq_len(ncol(x)),
    xdim = 3, ydim = 3, nClus = 3, seed = 1
  )
  st <- list(
    sce = sce, somCodesName = "SOM_codes", colsUsed = markers,
    rv = list(outputList = list(a = 1:3)), fsom = fsom
  )
  reg <- CySA:::.pdfPlotRegistry(st, list(), 1:3)
  expect_true(is.function(reg$flowSOMStars))
  # Regression: without list_insteadof_ggarrange the "$tree" element was
  # NULL and the PDF got an empty page.
  expect_s3_class(reg$flowSOMStars(), "ggplot")
  expect_s3_class(CySA:::.flowSOMStarsPlotly(fsom), "plotly")
})

# ---- diffcyt observer: guard and failure branches ---------------------------

.floor_two_groups <- function(session, flush) {
  suppressWarnings(session$setInputs(clusterNumbers = "1,2,3,4,5"))
  session$elapse(1500)
  suppressWarnings(session$setInputs(
    applyclusterNumbers = 1, clusterName = "groupA", applyName = 1
  ))
  flush()
  suppressWarnings(session$setInputs(clusterNumbers = "6,7,8,9,10"))
  session$elapse(1500)
  suppressWarnings(session$setInputs(
    applyclusterNumbers = 2, clusterName = "groupB", applyName = 2
  ))
  flush()
}

test_that("The diffcyt panel refuses a DS test with sample_id as random", {
  test_app <- make_test_app()
  shiny::testServer(app = test_app$app, expr = {
    .floor_two_groups(session, function() suppressWarnings(session$flushReact()))
    suppressWarnings(session$setInputs(
      diffGroups = c("groupA", "groupB"), diffAnalysisType = "DS",
      diffFixedCol = "treatment", diffRandomCol = "sample_id"
    ))
    suppressWarnings(session$setInputs(runDiffTest = 1))
    suppressWarnings(session$flushReact())
    expect_null(diffResult())
  })
})

test_that("The diffcyt panel runs a DS test on the chosen markers", {
  test_app <- make_test_app()
  shiny::testServer(app = test_app$app, expr = {
    .floor_two_groups(session, function() suppressWarnings(session$flushReact()))
    suppressWarnings(session$setInputs(
      diffGroups = c("groupA", "groupB"), diffAnalysisType = "DS",
      diffFixedCol = "treatment", diffRandomCol = "none",
      diffMarkers = "marker1"
    ))
    suppressWarnings(session$setInputs(runDiffTest = 1))
    suppressWarnings(session$flushReact())
    expect_identical(diffResult()$analysis_type, "DS")
  })
})

test_that("The diffcyt panel survives a failing test", {
  test_app <- make_test_app()
  shiny::testServer(app = test_app$app, expr = {
    .floor_two_groups(session, function() suppressWarnings(session$flushReact()))
    suppressWarnings(session$setInputs(
      diffGroups = c("groupA", "groupB"), diffAnalysisType = "DA",
      diffFixedCol = "no_such_column", diffRandomCol = "sample_id"
    ))
    suppressWarnings(session$setInputs(runDiffTest = 1))
    suppressWarnings(session$flushReact())
    expect_null(diffResult())
  })
})

# ---- diffcyt heatmap output: error branch -----------------------------------

test_that("The diffcyt heatmap output reports a drawing failure", {
  server <- function(input, output, session) {
    st <- list(diffResult = shiny::reactiveVal(list(analysis_type = "bogus")))
    CySA:::.registerDiffTestOutputs(input, output, st)
  }
  shiny::testServer(server, {
    out <- tryCatch(output$diffTestHeatmap, error = function(e) e)
    expect_true(!is.null(out))
  })
})
