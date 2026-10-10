test_that("clusterSelector() returns a shiny app object", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
  som_codes <- S4Vectors::metadata(sce)$SOM_codes

  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)),
    label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id
  )

  markers <- S4Vectors::metadata(sce)$map$colsUsed
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  }

  arr <- array(
    data = seq_len(10 * 10 * length(markers)),
    dim = c(10, 10, length(markers))
  )
  somRasterObj <- raster::brick(arr)
  names(somRasterObj) <- markers

  app <- clusterSelector(
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList,
    dend = dend,
    dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = somRasterObj
  )

  expect_s3_class(app, "shiny.appobj")
})

test_that("clusterSelector() validates missing SingleCellExperiment", {
  expect_error(
    clusterSelector(
      sce = "not an sce",
      sce_subsampled = "not an sce",
      dList = list(c("a", "b")),
      dend = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
      dendTable = data.frame(id = 1),
      clusterPatientTable = table(c(1)),
      somRasterData = data.frame(x = 1, y = 1)
    ),
    "must be a SingleCellExperiment"
  )
})

test_that("clusterSelector() validates required metadata", {
  sce <- CySA_example_sce()
  sce_bad <- sce
  S4Vectors::metadata(sce_bad)$map <- NULL

  expect_error(
    clusterSelector(
      sce = sce_bad,
      sce_subsampled = sce_bad,
      dList = list(c("marker1", "marker2")),
      dend = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
      dendTable = data.frame(id = 1),
      clusterPatientTable = table(c(1)),
      somRasterData = data.frame(x = 1, y = 1, marker1 = 1)
    ),
    "colsUsed is required"
  )
})

test_that("clusterSelector() initialises outputList in env", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
  som_codes <- S4Vectors::metadata(sce)$SOM_codes

  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)),
    label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id
  )

  markers <- S4Vectors::metadata(sce)$map$colsUsed
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  }

  env <- new.env()
  app <- clusterSelector(
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList,
    dend = dend,
    dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = NULL,
    env = env
  )

  expect_s3_class(app, "shiny.appobj")
  expect_true("outputList" %in% ls(envir = env))
  expect_true("Rest" %in% names(env$outputList))
})

test_that("clusterSelector() preserves user-supplied outputList groups", {
  sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
  som_codes <- S4Vectors::metadata(sce)$SOM_codes

  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)),
    label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id
  )

  markers <- S4Vectors::metadata(sce)$map$colsUsed
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) {
    somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  }

  env <- new.env()
  initial <- list(MyGroup = c(1L, 2L))
  app <- clusterSelector(
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled,
    outputList = initial,
    dList = prepped$dList,
    dend = dend,
    dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = NULL,
    env = env
  )

  expect_s3_class(app, "shiny.appobj")
  expect_true("MyGroup" %in% names(env$outputList))
  expect_true("Rest" %in% names(env$outputList))
  expect_true(all(env$outputList$MyGroup == c(1L, 2L)))
})

test_that(".buildFaviconLink returns a valid SVG data-URI <link> tag", {
  link <- CySA:::.buildFaviconLink()
  html <- htmltools::doRenderTags(link)

  expect_match(html, '<link rel="icon"', fixed = TRUE)
  expect_match(html, "data:image/svg\\+xml,")

  href <- sub('.*href="([^"]*)".*', "\\1", html)
  svg <- utils::URLdecode(sub("^data:image/svg\\+xml,", "", href))
  expect_match(svg, "^<svg xmlns=")
  # 4 blue nodes + 1 red "selected" node, the app's own selection colour.
  expect_equal(lengths(regmatches(svg, gregexpr("<circle", svg))), 5L)
})
