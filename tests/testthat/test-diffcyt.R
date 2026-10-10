# Tests for R/diffcyt-helpers.R ----
# diffcyt differential abundance (DA) / differential state (DS) testing on
# CySA's own saved cluster groups, adapted from the reference diffcyt
# workflow on CATALYST metaclusters reproduced in
# vignettes/catalyst_workflow.Rmd.

.diffTestFixture <- function(n_cells = 400, n_nodes = 20, n_markers = 3) {
  sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes,
    n_markers = n_markers)
  S4Vectors::metadata(sce)$experiment_info$treatment <- factor(
    rep(c("control", "treated"),
      length.out = nrow(S4Vectors::metadata(sce)$experiment_info)),
    levels = c("control", "treated")
  )
  outputList <- list(
    groupA = as.character(seq_len(n_nodes %/% 2)),
    groupB = as.character((n_nodes %/% 2 + 1):n_nodes)
  )
  list(sce = sce, outputList = outputList,
    experimentInfo = S4Vectors::metadata(sce)$experiment_info)
}

test_that(".nodeGroupMap maps every saved node id to its group name", {
  outputList <- list(groupA = c(1L, 2L, 3L), groupB = c("4", "5"))
  m <- CySA:::.nodeGroupMap(outputList)

  expect_equal(unname(m[c("1", "2", "3")]), rep("groupA", 3))
  expect_equal(unname(m[c("4", "5")]), rep("groupB", 2))
})

test_that(".checkGroupsDisjoint errors on overlapping groups", {
  expect_silent(CySA:::.checkGroupsDisjoint(
    list(groupA = c(1, 2), groupB = c(3, 4))
  ))
  expect_error(
    CySA:::.checkGroupsDisjoint(list(groupA = c(1, 2), groupB = c(2, 3))),
    "overlap"
  )
})

test_that("runDiffTest rejects overlapping groups instead of crashing", {
  # Regression test: "selected" (or any ad hoc group) sharing SOM nodes with
  # a named group used to make .nodeGroupMap() silently drop the second
  # group's cells, producing an all-NA row that crashed
  # diffcyt::plotHeatmap()'s internal hclust() row clustering.
  fx <- .diffTestFixture()
  overlapping <- list(
    groupA = fx$outputList$groupA, groupB = fx$outputList$groupA,
    Rest = fx$outputList$groupB
  )
  expect_error(
    runDiffTest(
      fx$sce, overlapping, fx$experimentInfo, fixedCol = "treatment",
      randomCol = "sample_id", analysisType = "DA"
    ),
    "overlap"
  )
})

test_that("runDiffTest drops an empty saved group with a warning", {
  fx <- .diffTestFixture()
  withGroups <- c(fx$outputList, list(emptyGroup = integer(0)))
  expect_warning(
    res <- runDiffTest(
      fx$sce, withGroups, fx$experimentInfo, fixedCol = "treatment",
      randomCol = "sample_id", analysisType = "DA"
    ),
    "Dropped empty"
  )
  expect_false("emptyGroup" %in%
    as.character(SummarizedExperiment::rowData(res$res)$cluster_id))
})

test_that("runDiffTest errors with fewer than 2 non-empty groups", {
  fx <- .diffTestFixture()
  expect_error(
    runDiffTest(
      fx$sce, list(onlyGroup = fx$outputList$groupA, emptyGroup = integer(0)),
      fx$experimentInfo, fixedCol = "treatment", randomCol = "sample_id",
      analysisType = "DA"
    ),
    "At least 2"
  )
})

test_that(".buildDiffcytSce drops cells outside every saved group", {
  fx <- .diffTestFixture()
  dSe <- CySA:::.buildDiffcytSce(
    fx$sce, fx$outputList, fx$experimentInfo,
    typeMarkers = rownames(fx$sce)
  )

  expect_true(all(c("sample_id", "cluster_id") %in%
    colnames(SummarizedExperiment::rowData(dSe))))
  expect_setequal(
    levels(SummarizedExperiment::rowData(dSe)$cluster_id),
    names(fx$outputList)
  )
  expect_equal(nrow(dSe), sum(fx$sce$cluster_id %in%
    as.integer(unlist(fx$outputList))))
})

test_that(".buildDiffcytSce errors when no cell belongs to a saved group", {
  fx <- .diffTestFixture()
  expect_error(
    CySA:::.buildDiffcytSce(
      fx$sce, list(empty = 999L), fx$experimentInfo,
      typeMarkers = rownames(fx$sce)
    ),
    "save at least"
  )
})

test_that("runDiffTest DA produces a non-NA p-value with a random effect", {
  fx <- .diffTestFixture()
  res <- runDiffTest(
    fx$sce, fx$outputList, fx$experimentInfo, fixedCol = "treatment",
    randomCol = "sample_id", analysisType = "DA"
  )

  expect_identical(res$analysis_type, "DA")
  tbl <- diffTestTable(res)
  expect_equal(nrow(tbl), length(fx$outputList))
  expect_false(anyNA(tbl$p_val))
})

test_that("runDiffTest DS produces results for the requested markers only", {
  fx <- .diffTestFixture()
  markers <- rownames(fx$sce)[1:2]
  res <- runDiffTest(
    fx$sce, fx$outputList, fx$experimentInfo, fixedCol = "treatment",
    randomCol = "none", analysisType = "DS", markers = markers
  )

  expect_identical(res$analysis_type, "DS")
  tbl <- diffTestTable(res)
  expect_setequal(as.character(unique(tbl$marker_id)), markers)
  expect_equal(nrow(tbl), length(fx$outputList) * length(markers))
})

test_that("plotDiffTestHeatmap draws a ComplexHeatmap for DA and DS", {
  fx <- .diffTestFixture()
  daRes <- runDiffTest(
    fx$sce, fx$outputList, fx$experimentInfo, fixedCol = "treatment",
    randomCol = "sample_id", analysisType = "DA"
  )
  dsRes <- runDiffTest(
    fx$sce, fx$outputList, fx$experimentInfo, fixedCol = "treatment",
    randomCol = "none", analysisType = "DS", markers = rownames(fx$sce)[1]
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s4_class(plotDiffTestHeatmap(daRes), "HeatmapList")
  expect_s4_class(plotDiffTestHeatmap(dsRes), "HeatmapList")
})

test_that("The diffcyt panel runs a DA test end to end inside the app", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    # clusterNumbers feeds the debounced inputClusterNumber() reactive
    # (1000ms, server-reactives.R) that applyName reads; each group needs
    # the debounce to elapse before the next one is typed, or applyName
    # saves a stale (previous) selection -- which, for two groups made back
    # to back, would make them accidentally share every node.
    suppressWarnings(session$setInputs(clusterNumbers = "1,2,3,4,5"))
    session$elapse(1500)
    suppressWarnings(session$setInputs(
      applyclusterNumbers = 1, clusterName = "groupA", applyName = 1
    ))
    suppressWarnings(session$flushReact())
    suppressWarnings(session$setInputs(clusterNumbers = "6,7,8,9,10"))
    session$elapse(1500)
    suppressWarnings(session$setInputs(
      applyclusterNumbers = 2, clusterName = "groupB", applyName = 2
    ))
    suppressWarnings(session$flushReact())
    expect_identical(rv$outputList$groupA, c(1L, 2L, 3L, 4L, 5L))
    expect_identical(rv$outputList$groupB, c(6L, 7L, 8L, 9L, 10L))

    suppressWarnings(session$setInputs(
      diffGroups = c("groupA", "groupB"), diffAnalysisType = "DA",
      diffFixedCol = "treatment", diffRandomCol = "sample_id"
    ))
    suppressWarnings(session$setInputs(runDiffTest = 1))
    suppressWarnings(session$flushReact())

    result <- diffResult()
    expect_false(is.null(result))
    expect_identical(result$analysis_type, "DA")
  })
})

test_that("The diffcyt panel refuses a DA test with no random effect", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1,2,3,4,5", applyclusterNumbers = 1,
      clusterName = "groupA", applyName = 1
    ))
    suppressWarnings(session$flushReact())

    suppressWarnings(session$setInputs(
      diffAnalysisType = "DA", diffFixedCol = "treatment",
      diffRandomCol = "none"
    ))
    suppressWarnings(session$setInputs(runDiffTest = 1))
    suppressWarnings(session$flushReact())

    expect_null(diffResult())
  })
})
