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

test_that("node_dispersion is appended once and is NA for empty nodes", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  isNode10 <- as.integer(sce$cluster_id) == 10L
  sce$cluster_id[isNode10] <- levels(sce$cluster_id)[1L]
  out <- .ensureSOMStatsDispersion(sce)
  st <- S4Vectors::metadata(out)$SOM_stats
  expect_true("node_dispersion" %in% colnames(st))
  expect_identical(nrow(st), 10L)
  expect_true(is.na(st$node_dispersion[10]))
  expect_true(all(is.finite(st$node_dispersion[1:9])))
  out2 <- .ensureSOMStatsDispersion(out)
  expect_identical(colnames(S4Vectors::metadata(out2)$SOM_stats), colnames(st))
})

test_that("node_dispersion attributes SDs by label value, not level order", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  levs <- as.character(seq_len(10))
  # Relevel the same labels: with a non-identity level order the level
  # INDEX of a cell's cluster_id no longer equals its label value, so a
  # level-index lookup would misattribute the per-node SDs.
  sce$cluster_id <-
    factor(levs[as.integer(sce$cluster_id)], levels = rev(levs))
  exprs <- SummarizedExperiment::assay(sce, "exprs")
  clusterId <- as.integer(as.character(
    SingleCellExperiment::colData(sce)$cluster_id))
  ref <- vapply(seq_len(10), function(lv) {
    idx <- which(clusterId == lv)
    if (length(idx) == 0L) {
      return(NA_real_)
    }
    mean(matrixStats::rowSds(exprs[, idx, drop = FALSE]))
  }, numeric(1))
  out <- .ensureSOMStatsDispersion(sce)
  got <- S4Vectors::metadata(out)$SOM_stats$node_dispersion
  expect_true(all(is.finite(got)))
  expect_equal(got, ref)
})

test_that("stale id-bearing cache labels are refreshed by id", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = rep(c("A", "B"), length.out = 10))
  somFile <- file.path(tempdir(),
    paste0("cysa-cache-", Sys.getpid(), "-stale-id.RDS"))
  cacheF <- sub("\\.RDS$", ".subsampled.RData", somFile)
  on.exit(unlink(c(somFile, cacheF)), add = TRUE)
  saveRDS(sce, somFile)
  sce2 <- sce
  sce2$id <- seq_len(ncol(sce2))
  # A cache from a run whose labels differ from the full object's current
  # ones: the id path must refresh them rather than trust the cache.
  sce_subsampled <- sce2[, 1:100]
  sce_subsampled$cell_type <- "STALE"
  save(sce_subsampled, file = cacheF)
  prepped <- prepClusterSelectorData(sce2, somFile = somFile,
    total_cells_to_sample = 100, cellTypeCol = "cell_type")
  sub <- prepped$sce_subsampled
  expect_false(all(as.character(
    SingleCellExperiment::colData(sub)$cell_type) == "STALE"))
  expect_identical(
    as.character(SingleCellExperiment::colData(sub)$cell_type),
    as.character(sce2$cell_type[
      SingleCellExperiment::colData(sub)$id]))
})

test_that("stale label column in a cache without id is redrawn", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10,
    cell_types = rep(c("A", "B"), length.out = 10))
  somFile <- file.path(tempdir(),
    paste0("cysa-cache-", Sys.getpid(), "-stale-noid.RDS"))
  cacheF <- sub("\\.RDS$", ".subsampled.RData", somFile)
  on.exit(unlink(c(somFile, cacheF)), add = TRUE)
  saveRDS(sce, somFile)
  # A cache carrying a stale label column but no id cannot be refreshed
  # in place, so the subsample is redrawn from the full object and the
  # labels ride along in its colData.
  sce_subsampled <- sce[, 1:100]
  sce_subsampled$cell_type <- "STALE"
  save(sce_subsampled, file = cacheF)
  prepped <- prepClusterSelectorData(sce, somFile = somFile,
    total_cells_to_sample = 100, cellTypeCol = "cell_type")
  sub <- prepped$sce_subsampled
  expect_false(all(as.character(
    SingleCellExperiment::colData(sub)$cell_type) == "STALE"))
  # The seeded redraw reproduces the same cells as the original draw.
  ref <- withr::with_seed(123, .subsampleSce(sce, 100, "SOM_codes"))
  expect_identical(SummarizedExperiment::colnames(sub),
    SummarizedExperiment::colnames(ref))
})
