test_that("prepClusterSelectorData() returns required components", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)

  expect_type(prepped, "list")
  expect_named(prepped, c("sce", "sce_subsampled", "dList"))
  expect_s4_class(prepped$sce, "SingleCellExperiment")
  expect_s4_class(prepped$sce_subsampled, "SingleCellExperiment")
  expect_type(prepped$dList, "list")
})

test_that("prepClusterSelectorData() subsampling reduces cell count", {
  sce <- CySA_example_sce(n_cells = 500, n_nodes = 10)
  prepped <- prepClusterSelectorData(
    sce,
    total_cells_to_sample = 100,
    seed = 42
  )
  expect_lt(ncol(prepped$sce_subsampled), ncol(sce))
  expect_true(ncol(prepped$sce_subsampled) >= 100)
})

test_that("prepClusterSelectorData() works with user-provided dList", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  dList <- list(
    d1 = c("marker1", "marker2"),
    d2 = c("marker3", "marker4")
  )
  prepped <- prepClusterSelectorData(
    sce,
    dList = dList,
    total_cells_to_sample = 100
  )
  expect_equal(prepped$dList, dList)
})

test_that("prepClusterSelectorData() is reproducible with seed", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  prepped1 <- prepClusterSelectorData(
    sce,
    total_cells_to_sample = 100,
    seed = 123
  )
  prepped2 <- prepClusterSelectorData(
    sce,
    total_cells_to_sample = 100,
    seed = 123
  )
  expect_equal(
    SingleCellExperiment::colData(prepped1$sce_subsampled)$id,
    SingleCellExperiment::colData(prepped2$sce_subsampled)$id
  )
})

test_that("prepClusterSelectorData() rejects invalid total_cells_to_sample", {
  sce <- CySA_example_sce(n_cells = 100, n_nodes = 10)
  expect_error(
    prepClusterSelectorData(sce, total_cells_to_sample = -1),
    "must be a positive scalar"
  )
  expect_error(
    prepClusterSelectorData(sce, total_cells_to_sample = "many"),
    "must be a positive scalar"
  )
})

test_that("cellTypeCol adds node_labels with purity and majority", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("type", rep(1:5, length.out = 10)))
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  nl <- S4Vectors::metadata(prepped$sce)$node_labels
  expect_s3_class(nl, "data.frame")
  expect_identical(colnames(nl), c("id", "purity", "cell_type_majority"))
  expect_identical(nrow(nl), 10L)
  expect_true(all(nl$purity == 1))
  expect_setequal(na.omit(nl$cell_type_majority), paste0("type", 1:5))
})

test_that("NA labels: purity counts labeled cells only; empty nodes stay NA", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = paste0("type", rep(1:5, length.out = 10)))
  lab <- sce$cell_type
  lab[sce$cluster_id %in% c(9, 10)] <- NA
  sce$cell_type <- lab
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type")
  nl <- S4Vectors::metadata(prepped$sce)$node_labels
  expect_true(all(is.na(nl$purity[9:10])))
  expect_true(all(nl$purity[1:8] == 1))
})

test_that("unknown cellTypeCol stops", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  expect_error(prepClusterSelectorData(sce, total_cells_to_sample = 100,
    cellTypeCol = "nope"), "cellTypeCol")
})

test_that("cached subsample without the label column gets labels merged", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = rep(c("A", "B"), length.out = 10))
  somFile <- file.path(tempdir(),
    paste0("cysa-cache-", Sys.getpid(), ".RDS"))
  cacheF <- sub("\\.RDS$", ".subsampled.RData", somFile)
  on.exit(unlink(c(somFile, cacheF)), add = TRUE)
  saveRDS(sce, somFile)
  # A cache written before the label column existed: no cell_type in colData.
  # (save() must be given the object by name: "save(x = y)" saves nothing.)
  sce_subsampled <- sce[, 1:100]
  save(sce_subsampled, file = cacheF)
  prepped <- prepClusterSelectorData(sce, somFile = somFile,
    total_cells_to_sample = 100, cellTypeCol = "cell_type")
  sub <- prepped$sce_subsampled
  expect_true("cell_type" %in% colnames(SingleCellExperiment::colData(sub)))
  # With an id column the labels merge in place instead of re-subsampling.
  sce2 <- sce
  sce2$id <- seq_len(ncol(sce2))
  sce_subsampled <- sce2[, 1:100]
  save(sce_subsampled, file = cacheF)
  prepped2 <- prepClusterSelectorData(sce2, somFile = somFile,
    total_cells_to_sample = 100, cellTypeCol = "cell_type")
  env <- new.env()
  load(cacheF, envir = env)
  expect_identical(
    as.character(SingleCellExperiment::colData(prepped2$sce_subsampled)$cell_type),
    as.character(sce2$cell_type[env$sce_subsampled$id]))
})
