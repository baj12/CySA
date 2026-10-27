test_that("CySA_example_sce() returns a valid SingleCellExperiment", {
  sce <- CySA_example_sce()

  expect_s4_class(sce, "SingleCellExperiment")
  expect_true("exprs" %in% names(SummarizedExperiment::assays(sce)))
  expect_true("cluster_id" %in% names(SingleCellExperiment::colData(sce)))
  expect_true("sample_id" %in% names(SingleCellExperiment::colData(sce)))

  metaD <- S4Vectors::metadata(sce)
  expect_true("SOM_codes" %in% names(metaD))
  expect_true("SOM_stats" %in% names(metaD))
  expect_true("map" %in% names(metaD))
  expect_true("colsUsed" %in% names(metaD$map))

  expect_equal(nrow(SummarizedExperiment::assay(sce, "exprs")), 12)
  expect_equal(ncol(sce), 1000)
})

test_that("CySA_example_sce() respects n_cells and n_nodes", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10, n_markers = 8)
  expect_equal(ncol(sce), 200)
  expect_equal(nrow(S4Vectors::metadata(sce)$SOM_codes), 10)
  expect_equal(nrow(SummarizedExperiment::assay(sce, "exprs")), 8)
})

test_that("CySA_example_sce adds a cell_type column from cell_types", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("type", rep(1:5, length.out = 10)))
  expect_true("cell_type" %in% colnames(SingleCellExperiment::colData(sce)))
  cd <- as.data.frame(SingleCellExperiment::colData(sce))
  expect_identical(as.character(cd$cell_type),
    paste0("type", rep(1:5, length.out = 10))[as.integer(cd$cluster_id)])
})

test_that("CySA_example_sce without cell_types is unchanged", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  expect_false("cell_type" %in% colnames(SingleCellExperiment::colData(sce)))
})

test_that("CySA_default_cluster_cols() returns a non-empty character vector", {
  cols <- CySA_default_cluster_cols()
  expect_type(cols, "character")
  expect_true(length(cols) >= 1)
  expect_true(all(grepl("^#", cols)))
})
