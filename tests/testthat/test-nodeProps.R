# Tests for node-level properties: nodeAmbiguity(), plotNodeProps() ----

test_that("nodeAmbiguity returns one row per populated node", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  amb <- nodeAmbiguity(sce, maxCells = 100, seed = 1)
  expect_s3_class(amb, "data.frame")
  expect_identical(colnames(amb),
    c("id", "n", "meanDistToCode", "medianRatio", "pctAmbiguous"))
  expect_setequal(amb$id, seq_len(10))
  # The squared-distance cancellation can push the ratio a hair above 1
  # (same artifact as computeClusterAmbiguity); allow a small tolerance.
  expect_true(all(amb$medianRatio > -1e-9 & amb$medianRatio <= 1 + 1e-3))
  expect_true(all(amb$pctAmbiguous >= 0 & amb$pctAmbiguous <= 100))
})

test_that("nodeAmbiguity caps cells and is seed-stable", {
  sce <- CySA_example_sce(n_cells = 500, n_nodes = 10)
  a <- nodeAmbiguity(sce, maxCells = 60, seed = 1)
  b <- nodeAmbiguity(sce, maxCells = 60, seed = 1)
  expect_identical(a, b)
  expect_true(sum(a$n) <= 60)
})

test_that("computeClusterAmbiguity still works after the refactor", {
  sce <- CySA_example_sce(n_cells = 100, n_nodes = 10)
  s <- computeClusterAmbiguity(sce, S4Vectors::metadata(sce)$SOM_codes)
  expect_identical(colnames(s),
    c("cluster_id", "n", "meanDistToCode", "medianRatio", "pctAmbiguous"))
})

test_that("plotNodeProps defaults to purity x n when labels exist", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("t", rep(1:2, 5)))
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  amb <- nodeAmbiguity(prepped$sce_subsampled, maxCells = 100)
  p <- plotNodeProps(prepped$sce, ambiguity = amb)
  expect_s3_class(p, "ggplot")
  expect_identical(as.character(p$labels$x), "purity")
  expect_identical(as.character(p$labels$y), "n")
})

test_that("plotNodeProps defaults to medianRatio x n without labels", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  p <- plotNodeProps(sce, ambiguity = nodeAmbiguity(sce, maxCells = 100))
  expect_identical(as.character(p$labels$x), "medianRatio")
  expect_identical(as.character(p$labels$y), "n")
})

test_that("plotNodeProps stops with too few numeric properties", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  amb <- nodeAmbiguity(sce, maxCells = 100)
  expect_error(plotNodeProps(sce, xVar = "notThere", ambiguity = amb),
    "notThere")
})
