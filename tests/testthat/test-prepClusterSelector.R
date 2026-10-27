# Tests for R/prepClusterSelector.R helpers and main function ----

# ══════════════════════════════════════════════════════════════════════════════
# .defaultDList
# ══════════════════════════════════════════════════════════════════════════════

test_that(".defaultDList builds six pairs from colsUsed", {
  sce <- CySA_example_sce(n_markers = 12)
  dl <- CySA:::.defaultDList(sce)
  expect_type(dl, "list")
  expect_length(dl, 6L)
  expect_true(all(vapply(dl, length, integer(1)) == 2L))
})

test_that(".defaultDList falls back to rownames when colsUsed too short", {
  sce <- CySA_example_sce(n_markers = 12)
  S4Vectors::metadata(sce)$map$colsUsed <- S4Vectors::metadata(sce)$map$colsUsed[1:10]
  dl <- CySA:::.defaultDList(sce)
  expect_equal(length(dl), 6L)
})

test_that(".defaultDList errors when fewer than 12 names available", {
  sce <- CySA_example_sce(n_markers = 6)
  S4Vectors::metadata(sce)$map$colsUsed <- NULL
  expect_error(CySA:::.defaultDList(sce), "at least 12")
})

# ══════════════════════════════════════════════════════════════════════════════
# prepClusterSelectorData
# ══════════════════════════════════════════════════════════════════════════════

test_that("prepClusterSelectorData uses cache file when present", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 6)
  sce_subsampled <- sce[, 1:50]
  tmp <- tempfile(fileext = ".RData")
  save(sce_subsampled = sce_subsampled, file = tmp)
  cache <- sub("\\.RData$", ".subsampled.RData", tmp)
  file.copy(tmp, cache, overwrite = TRUE)

  prepped <- prepClusterSelectorData(sce, somFile = tmp, total_cells_to_sample = 50)
  expect_s4_class(prepped$sce_subsampled, "SingleCellExperiment")
  expect_equal(ncol(prepped$sce_subsampled), 50L)

  unlink(c(tmp, cache))
})

# ══════════════════════════════════════════════════════════════════════════════
# .buildSOMCodes
# ══════════════════════════════════════════════════════════════════════════════

test_that(".buildSOMCodes recomputes missing markers", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 6, n_markers = 8)
  # Drop two markers from SOM_codes so they have to be recomputed.
  keep <- colnames(S4Vectors::metadata(sce)$SOM_codes)[1:6]
  S4Vectors::metadata(sce)$SOM_codes <- S4Vectors::metadata(sce)$SOM_codes[, keep, drop = FALSE]

  sce2 <- CySA:::.buildSOMCodes(sce)
  expect_true(all(rownames(sce2) %in% colnames(S4Vectors::metadata(sce2)$SOM_codes)))
})

test_that(".buildSOMCodes returns early when all markers present", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4)
  sce2 <- CySA:::.buildSOMCodes(sce)
  expect_equal(
    colnames(S4Vectors::metadata(sce)$SOM_codes),
    colnames(S4Vectors::metadata(sce2)$SOM_codes)
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# .buildSOMStats
# ══════════════════════════════════════════════════════════════════════════════

test_that(".buildSOMStats returns a data.frame with required columns", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 6)
  st <- CySA:::.buildSOMStats(sce)
  expect_s3_class(st, "data.frame")
  expect_true(all(c("median", "mean", "rdQu", "max", "n", "id") %in% colnames(st)))
  expect_equal(nrow(st), 6L)
})

test_that(".buildSOMStats handles empty clusters", {
  sce <- CySA_example_sce(n_cells = 100, n_nodes = 10)
  # force some clusters to have no cells
  cd <- SingleCellExperiment::colData(sce)
  cd$cluster_id <- factor(rep(1:2, length.out = ncol(sce)), levels = seq_len(10))
  SummarizedExperiment::colData(sce) <- cd
  st <- CySA:::.buildSOMStats(sce)
  expect_equal(nrow(st), 10L)
})

# ══════════════════════════════════════════════════════════════════════════════
# .buildSOMClusterStats
# ══════════════════════════════════════════════════════════════════════════════

test_that(".buildSOMClusterStats returns long data frame", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 4, n_markers = 6)
  df <- CySA:::.buildSOMClusterStats(sce, verbose = FALSE)
  expect_true(all(c("cluster", "stat", "variable", "value") %in% colnames(df)))
  expect_equal(length(unique(df$stat)), 6L)
})

# ══════════════════════════════════════════════════════════════════════════════
# .buildSOMRasterData
# ══════════════════════════════════════════════════════════════════════════════

test_that(".buildSOMRasterData returns square grid data", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 16, n_markers = 6)
  rd <- CySA:::.buildSOMRasterData(sce)
  expect_s3_class(rd, "data.frame")
  expect_true(all(c("x", "y", "N") %in% colnames(rd)))
})

test_that(".buildSOMRasterData errors for non-square grid", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10, n_markers = 6)
  expect_error(CySA:::.buildSOMRasterData(sce), "not square")
})

# ══════════════════════════════════════════════════════════════════════════════
# .buildDendTable
# ══════════════════════════════════════════════════════════════════════════════

test_that(".buildDendTable returns list without tree when useColTree = FALSE", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 6)
  res <- CySA:::.buildDendTable(sce, useColTree = FALSE)
  expect_type(res, "list")
  expect_named(res, c("dend", "dendTable", "colTree"))
  expect_null(res$colTree)
})

# ══════════════════════════════════════════════════════════════════════════════
# prepClusterSelector main function
# ══════════════════════════════════════════════════════════════════════════════

test_that("prepClusterSelector creates sce.shiny.RData", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 16, n_markers = 12)
  root <- tempfile("cysa")
  dir.create(file.path(root, "_data"), recursive = TRUE)
  saveRDS(sce, file.path(root, "_data", "sce.RDS"))

  res <- prepClusterSelector(root, total_cells_to_sample = 100, force = TRUE, verbose = FALSE)
  expect_equal(res, root)
  expect_true(file.exists(file.path(root, "_data", "sce.shiny.RData")))

  unlink(root, recursive = TRUE)
})

test_that("prepClusterSelector reuses existing output unless force = TRUE", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 16)
  root <- tempfile("cysa")
  dir.create(file.path(root, "_data"), recursive = TRUE)
  saveRDS(sce, file.path(root, "_data", "sce.RDS"))

  prepClusterSelector(root, total_cells_to_sample = 100, force = TRUE, verbose = FALSE)
  mtime1 <- file.info(file.path(root, "_data", "sce.shiny.RData"))$mtime

  Sys.sleep(1.1)
  prepClusterSelector(root, force = FALSE, verbose = FALSE)
  mtime2 <- file.info(file.path(root, "_data", "sce.shiny.RData"))$mtime

  expect_equal(mtime1, mtime2)
  unlink(root, recursive = TRUE)
})

test_that("prepClusterSelector errors when sce.RDS missing", {
  root <- tempfile("cysa")
  dir.create(file.path(root, "_data"), recursive = TRUE)
  expect_error(prepClusterSelector(root), "SCE file not found")
  unlink(root, recursive = TRUE)
})

test_that("prepClusterSelector uses provided subsampled SCE file", {
  sce <- CySA_example_sce(n_cells = 300, n_nodes = 16, n_markers = 12)
  root <- tempfile("cysa")
  dir.create(file.path(root, "_data"), recursive = TRUE)
  saveRDS(sce, file.path(root, "_data", "sce.RDS"))

  sub_file <- file.path(root, "sce.subsampled.RDS")
  saveRDS(sce[, 1:100], sub_file)

  prepClusterSelector(root, subsampledSCEFile = sub_file, force = TRUE, verbose = FALSE)
  expect_true(file.exists(file.path(root, "_data", "sce.shiny.RData")))
  unlink(root, recursive = TRUE)
})

test_that("prepClusterSelector bundle gains node_dispersion (build and reuse)", {
  # Build branch: no SOM_stats in the source SCE.
  sce1 <- CySA_example_sce(n_cells = 300, n_nodes = 16)
  S4Vectors::metadata(sce1)$SOM_stats <- NULL
  root1 <- tempfile("cysa")
  dir.create(file.path(root1, "_data"), recursive = TRUE)
  saveRDS(sce1, file.path(root1, "_data", "sce.RDS"))
  prepClusterSelector(root1, total_cells_to_sample = 100, force = TRUE,
    verbose = FALSE)
  env1 <- new.env()
  load(file.path(root1, "_data", "sce.shiny.RData"), envir = env1)
  expect_true("node_dispersion" %in%
    colnames(S4Vectors::metadata(env1$sce)$SOM_stats))
  unlink(root1, recursive = TRUE)

  # Reuse branch: existing SOM_stats (the example's, without the column)
  # gains it via prepClusterSelector.
  sce2 <- CySA_example_sce(n_cells = 300, n_nodes = 16)
  expect_false("node_dispersion" %in%
    colnames(S4Vectors::metadata(sce2)$SOM_stats))
  root2 <- tempfile("cysa")
  dir.create(file.path(root2, "_data"), recursive = TRUE)
  saveRDS(sce2, file.path(root2, "_data", "sce.RDS"))
  prepClusterSelector(root2, total_cells_to_sample = 100, force = TRUE,
    verbose = FALSE)
  env2 <- new.env()
  load(file.path(root2, "_data", "sce.shiny.RData"), envir = env2)
  expect_true("node_dispersion" %in%
    colnames(S4Vectors::metadata(env2$sce)$SOM_stats))
  unlink(root2, recursive = TRUE)
})

test_that("prepClusterSelector forwards cellTypeCol on the normal path", {
  root <- tempfile("cysa")
  dir.create(file.path(root, "_data"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 16,
    cell_types = paste0("t", rep(1:4, length.out = 16)))
  saveRDS(sce, file.path(root, "_data", "sce.RDS"))
  prepClusterSelector(root, sce = sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type", verbose = FALSE)
  env <- new.env()
  load(file.path(root, "_data", "sce.shiny.RData"), envir = env)
  expect_true(!is.null(S4Vectors::metadata(env$sce)$node_labels))
  expect_true("node_dispersion" %in%
    colnames(S4Vectors::metadata(env$sce)$SOM_stats))
  expect_true("cell_type" %in%
    colnames(SummarizedExperiment::colData(env$sce_subsampled)))
})

test_that("early return still applies label stats to the existing bundle", {
  root <- tempfile("cysa")
  dir.create(file.path(root, "_data"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 16,
    cell_types = paste0("t", rep(1:4, length.out = 16)))
  saveRDS(sce, file.path(root, "_data", "sce.RDS"))
  # First run WITHOUT labels (creates the bundle), then re-run WITH
  # cellTypeCol and force = FALSE: the early-return path must apply the
  # label stats to the loaded bundle anyway.
  prepClusterSelector(root, sce = sce, total_cells_to_sample = 100,
    verbose = FALSE)
  prepClusterSelector(root, sce = sce, total_cells_to_sample = 100,
    cellTypeCol = "cell_type", verbose = FALSE)
  env <- new.env()
  load(file.path(root, "_data", "sce.shiny.RData"), envir = env)
  expect_true(!is.null(S4Vectors::metadata(env$sce)$node_labels))
  expect_true("cell_type" %in%
    colnames(SummarizedExperiment::colData(env$sce_subsampled)))
})
