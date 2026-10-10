# tests/testthat/test-outputs_clusterSelector-gaps.R
#
# Targets branches uncovered ("!") in the coverage report for
# R/outputs_clusterSelector.R, driven through the real app from
# make_test_app() (helper-*.R) plus a local variant that adds a factor
# grouping column, which the stock CySA_example_sce() fixture lacks.
#
# Selection wiring (see observers_clusterSelector.R / server_clusterSelector.R):
#   input$clusterNumbers (text, comma-separated ids)
#     -> inputClusterNumber() [debounce 1000ms]
#     -> observer sets rsUsed(valid)
#     -> rsUsed_d = rsUsed %>% debounce(1500ms)
#
# rsUsed() feeds: scatter, dendPlot, somRasterPlot, cellCounts,
#   cellPercentages, countBarPlot, PercentBarPlot, flowSOMPiePlot.
# rsUsed_d() feeds: somData<i>, somDataMain, tsnePlot, umapPlot, pcaPlot.
#
# dimSelection() is populated from input$currentDimX/currentDimY via a
# 250ms-debounced activeDims() reactive, and is ALWAYS length 1 regardless
# of nPlots (see note in chat about the downloadPlots loop).

library(testthat)
library(shiny)

# helper-*.R (make_test_app) is auto-sourced by testthat before this file.

# ---- Local fixture: adds a usable factor grouping column -------------------
# CySA_example_sce()'s experiment_info only has sample_id (character) and
# some_numeric (numeric) -- no factor column, so groupsVar/group1/group2 and
# the real t.test() path are unreachable with the stock fixture.

make_test_app_with_groups <- function(n_cells = 400, n_nodes = 10) {
  set.seed(12)
  sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes)

  ei <- S4Vectors::metadata(sce)$experiment_info
  # sample1/sample2 -> "A", sample3/sample4 -> "B", each with >1 sample so
  # sd(x) and sd(y) aren't both zero in ttestResult.
  ei$cohort <- factor(rep(c("A", "A", "B", "B"), length.out = nrow(ei)))
  S4Vectors::metadata(sce)$experiment_info <- ei

  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
  som_codes <- S4Vectors::metadata(sce)$SOM_codes
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)), label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id, cluster_id = sce$cluster_id
  )
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  arr <- array(seq_len(10 * 10 * length(markers)), dim = c(10, 10, length(markers)))
  somRasterObj <- raster::brick(arr)
  names(somRasterObj) <- markers

  app <- clusterSelector(
    sce = prepped$sce, sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList, dend = dend, dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData, somRasterObj = somRasterObj,
    env = new.env()
  )
  list(app = app, sce = prepped$sce, markers = markers)
}

# Common: get past clusterNumbers -> rsUsed -> rsUsed_d chain.
select_nodes <- function(session, ids = "1, 2, 3", wait_for_debounced = TRUE) {
  session$setInputs(clusterNumbers = ids)
  session$elapse(1001) # inputClusterNumber() debounce
  if (wait_for_debounced) session$elapse(1501) # rsUsed_d debounce
}

# Common: get past dimSelection()'s req via currentDimX/currentDimY.
select_dims <- function(session, x, y) {
  session$setInputs(currentDimX = x, currentDimY = y)
  session$elapse(300) # activeDims() debounce (250ms)
}

# ---- Gap 1: scatter body past req() (real quantile/plotly logic) ----------

test_that("scatter renders a real plot once selection/dims/samples are valid", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_dims(session, "marker1", "marker2")
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(samples2plot = c("sample1", "sample2", "sample3", "sample4"))
      out <- output$scatter
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 2: somData1 -- both showGroups branches ---------------------------

test_that("somData1 uses the non-grouped path when showGroups is FALSE", {
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    select_dims(session, "marker1", "marker2")
    select_nodes(session, "1, 2, 3")
    session$setInputs(showGroups = FALSE, colorbyGroups = character(0))
    out <- output$somData1
    expect_true(!is.null(out))
  }))
})

test_that("somData1 uses the projection path when showGroups is TRUE", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_dims(session, "marker1", "marker2")
      select_nodes(session, "1, 2, 3")
      session$setInputs(showGroups = TRUE, colorbyGroups = "selected")
      out <- output$somData1
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 3: somDataMain ------------------------------------------------------

test_that("SOM 2D box renders the optional neighbourhood-agreement checkbox", {
  # The agreement view is opt-in: the checkbox exists, default unchecked.
  # shiny.appobj keeps the ui tag in the closure of serverFuncSource.
  app <- make_test_app()
  e <- environment(app$app$serverFuncSource)
  html <- as.character(e$ui)
  expect_true(grepl('id="nnAgreement"', html, fixed = TRUE))
  expect_true(grepl("Neighbourhood agreement", html, fixed = TRUE))
  # The checkbox input itself carries no checked attribute (later unrelated
  # checkboxes do, so anchor on the input tag).
  seg <- substr(html, regexpr('id="nnAgreement"', html, fixed = TRUE),
    nchar(html))
  nnInput <- regmatches(seg, regexpr("^id=\"nnAgreement\"[^>]*", seg))[[1L]]
  expect_false(grepl("checked", nnInput))
})

test_that("somDataMain renders past the dims/rsUsed_d guard", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_dims(session, "marker1", "marker2")
      select_nodes(session, "1, 2, 3")
      session$setInputs(showGroups = FALSE)
      out <- output$somDataMain
      expect_true(!is.null(out))
    })
  )
})

test_that("somDataMain carries agreement colours only when checked", {
  # The agreement colour scale must appear in the somDataMain plotly only
  # when nnAgreement is checked; group colouring or the plain scale must
  # never carry it.
  fx <- make_test_app()
  has_agreement_colourbar <- function(q) {
    if (is.character(q)) {
      q <- jsonlite::fromJSON(q, simplifyVector = FALSE)
    }
    any(vapply(q$x$data, function(tr) {
      cb <- tr$marker$colorbar %||% NULL
      !is.null(cb) && identical(cb$title, "agreement")
    }, logical(1)))
  }
  check <- function(nn, expected, showGroups = TRUE) {
    quiet_plotly_test(
      testServer(fx$app, {
        select_dims(session, "marker1", "marker2")
        select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
        session$setInputs(showGroups = showGroups, nnAgreement = nn,
          colorbyGroups = character(0))
        session$flushReact()
        session$elapse(1100) # somBasePlot's bindCache debounce
        out <- output$somDataMain
        expect_true(!is.null(out))
        expect_identical(has_agreement_colourbar(out), expected)
      })
    )
  }
  expect_error(check(TRUE, TRUE, showGroups = FALSE), regexp = NA)
  expect_error(check(TRUE, TRUE, showGroups = TRUE), regexp = NA)
  expect_error(check(FALSE, FALSE, showGroups = FALSE), regexp = NA)
})

# ---- Gap 4: tsne/umap/pca -- both showlegend branches, real computation ----

test_that("tsne/umap/pca compute and render with legend shown and hidden", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3")
      session$setInputs(
        dimRedSelection = fx$sce %>% (function(s) S4Vectors::metadata(s)$map$colsUsed)(),
        perplexity = 5, n_neighbors = 5, showlegend = TRUE
      )
      session$elapse(2100) # dimRedSelection debounce (1000) + tsne/umap/pca debounce (1000)

      expect_true(!is.null(output$tsne))
      expect_true(!is.null(output$umap))
      expect_true(!is.null(output$pca))

      session$setInputs(showlegend = FALSE)
      expect_true(!is.null(output$tsne))
      expect_true(!is.null(output$umap))
      expect_true(!is.null(output$pca))
    })
  )
})

# ---- Gap 5: cellCounts -- both comparison branches -------------------------

test_that("cellCounts adds a comparison row when compareStatsTo matches experiment_info", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(compareStatsTo = "some_numeric", singleNode = 1)
      tbl <- output$cellCounts
      expect_true(!is.null(tbl))
    })
  )
})

test_that("cellCounts adds a comparison row when compareStatsTo matches an outputList name", {
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
    session$setInputs(clusterName = "popA")
    session$setInputs(applyName = 1) # actionButton click
    session$setInputs(compareStatsTo = "popA", singleNode = 1)
    tbl <- output$cellCounts
    expect_true(!is.null(tbl))
  }))
})

test_that("cellCounts and cellPercentages survive a stale node id typed into a named group", {
  # Regression: "Apply name" saved the raw typed node ids (including any
  # out-of-range id) into outputList without validation; every output that
  # later indexed clusterPatientTable/SOM codes by that group then crashed
  # with "subscript out of bounds" (PercentBar, CountBar, cellCounts,
  # cellPercentages, VlnPlot, VlnPlot2 were all affected).
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    select_nodes(session, "1, 2, 99999", wait_for_debounced = FALSE)
    session$setInputs(clusterName = "grpBad")
    session$setInputs(applyName = 1)
    session$flushReact()
    expect_false(99999L %in% rv$outputList[["grpBad"]])

    session$setInputs(compareStatsTo = "grpBad", relativeTo = "grpBad",
      singleNode = 1,
      violinSelection = S4Vectors::metadata(fx$sce)$map$colsUsed)
    session$elapse(1200) # violinPlotSelection() debounce
    session$flushReact()
    expect_no_error(force(output$cellCounts))
    expect_no_error(force(output$cellPercentages))
    expect_no_error(force(output$CountBar))
    expect_no_error(force(output$PercentBar))
    expect_no_error(force(output$VlnPlot))
    expect_no_error(force(output$VlnPlot2))
  }))
})

test_that("cellPercentages works when relativeTo is a numeric experiment_info column", {
  # Regression: the inline eI matrix built by .cellPercentages() lost its
  # sample_id row names (unlike the shared .numericExperimentInfo()
  # helper), so indexing it by sample id in .computeRelativeCounts() threw
  # "subscript out of bounds".
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
    session$setInputs(relativeTo = "some_numeric")
    expect_no_error(force(output$cellPercentages))
  }))
})

# ---- Gap 6: ttestResult -- empty-data guard AND real t.test() path ---------

test_that("ttestResult reports 'no data' when groupsVar is 'none'", {
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    session$setInputs(groupsVar = "none", relativeTo = "none")
    select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
    expect_output(print(output$ttestResult), "no data")
  }))
})

test_that("outputs: ttestResult runs a real t-test when a valid grouping is set", {
    set.seed(42)

    n_markers   <- 6
    n_som_nodes <- 10
    n_samples   <- 4
    n_cells     <- 60  # per sample

    marker_names <- paste0("marker", seq_len(n_markers))
    sample_ids   <- paste0("S", seq_len(n_samples))
    conditions   <- factor(
        c("control", "control", "treated", "treated"),
        levels = c("control", "treated")
    )

    experiment_info <- data.frame(
        sample_id   = sample_ids,
        treatment   = conditions,
        total_cells = rep(1000L, n_samples),
        stringsAsFactors = FALSE
    )

    som_codes <- matrix(
        runif(n_som_nodes * n_markers),
        nrow = n_som_nodes, ncol = n_markers,
        dimnames = list(paste0("node", seq_len(n_som_nodes)), marker_names)
    )

    # Control samples draw mostly from nodes 1-5; treated mostly from 6-10,
    # with independent random draws per sample so per-sample counts differ.
    assign_clusters <- function(is_control, n) {
        pool <- if (is_control) 1:5 else 6:10
        c(
            base::sample(pool, round(n * 0.7), replace = TRUE),
            base::sample(seq_len(n_som_nodes), n - round(n * 0.7), replace = TRUE)
        )
    }

    coldata_list <- base::mapply(function(sid, cond) {
        data.frame(
            sample_id  = sid,
            cluster_id = assign_clusters(cond == "control", n_cells),
            stringsAsFactors = FALSE
        )
    }, sample_ids, as.character(conditions), SIMPLIFY = FALSE)
    coldata_df <- do.call(rbind, coldata_list)
    n_total <- nrow(coldata_df)

    exprs_mat <- matrix(
        runif(n_markers * n_total),
        nrow = n_markers,
        dimnames = list(marker_names, paste0("cell", seq_len(n_total)))
    )
    node_counts <- tabulate(coldata_df$cluster_id, nbins = n_som_nodes)

    som_stats <- data.frame(
        id     = seq_len(n_som_nodes),
        n      = node_counts,
        mean   = seq_len(n_som_nodes) / n_som_nodes,
        median = seq_len(n_som_nodes) / n_som_nodes,
        rdQu   = seq_len(n_som_nodes) / n_som_nodes,
        max    = seq_len(n_som_nodes) / n_som_nodes,
        stringsAsFactors = FALSE
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(exprs = exprs_mat),
        colData = S4Vectors::DataFrame(
            sample_id  = coldata_df$sample_id,
            cluster_id = factor(coldata_df$cluster_id, levels = seq_len(n_som_nodes))
        ),
        metadata = list(
            SOM_codes       = som_codes,
            SOM_stats       = som_stats,
            map             = list(colsUsed = marker_names),
            experiment_info = experiment_info
        )
    )
    rownames(sce) <- marker_names

    clusterPatientTable <- table(
        sample_id  = sce$sample_id,
        cluster_id = sce$cluster_id
    )

    dList <- list(
        d1 = c("marker1", "marker2"),
        d2 = c("marker3", "marker4"),
        d3 = c("marker5", "marker6"),
        d4 = c("marker1", "marker3"),
        d5 = c("marker2", "marker5"),
        d6 = c("marker4", "marker6")
    )

    prepped <- prepClusterSelectorData(
        sce,
        dList = dList,
        total_cells_to_sample = 200
    )
    dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
    dendTable <- data.frame(
        id = seq_len(n_som_nodes),
        label = rownames(som_codes),
        stringsAsFactors = FALSE
    )
    somRasterData <- data.frame(
        x  = rep(seq_len(5), length.out = n_som_nodes),
        y  = rep(seq_len(2), each = ceiling(n_som_nodes / 2)),
        id = seq_len(n_som_nodes)
    )
    for (m in marker_names) {
        somRasterData[[m]] <- som_codes[, m]
    }
    app <- clusterSelector(
        sce                 = prepped$sce,
        sce_subsampled      = prepped$sce_subsampled,
        dList               = prepped$dList,
        dend                = dend,
        dendTable           = dendTable,
        clusterPatientTable = clusterPatientTable,
        somRasterData       = somRasterData,
        somRasterObj        = NULL
    )

    shiny::testServer(
        app = app,
        expr = {
            suppressWarnings(session$setInputs(clusterNumbers = "1,2,3,4,5"))
            session$elapse(1000)
            suppressWarnings(session$flushReact())

            suppressWarnings(session$setInputs(
                groupsVar  = "treatment",
                group1     = "control",
                group2     = "treated",
                relativeTo = "none"
            ))
            session$elapse(1000)
            suppressWarnings(session$flushReact())

            result <- paste(output$ttestResult, collapse = "\n")

            expect_false(grepl("no data", result, fixed = TRUE))
            expect_false(grepl("not enough data", result, fixed = TRUE))
            expect_true(grepl("p-value", result, fixed = TRUE))
            expect_true(grepl("t = ", result, fixed = TRUE))
        }
    )
})


test_that("drawProjection produces non-empty traces for a real selection", {
  fx <- make_test_app()
  sce <- fx$sce
  pp1 <- plotSOMScatter(
    x = sce, chs = c("marker1", "marker2"),
    pointSize = "max", color_by = "n", zeros = TRUE,
    xRN = rownames(sce), xCN = colnames(sce)
  )
  dimSel <- list(list(dims = c("marker1", "marker2")))
  projectionDf <- buildProjectionDf(pp1, 1L, dimSel, sce, somCodesName = "SOM_codes")

  p <- drawProjection(projectionDf,
    rs = 1:3, colorbyGroups = "selected",
    sce = sce, outputList = list(selected = 1:3)
  )

  built <- plotly::plotly_build(plotly::ggplotly(p))
  trace_lengths <- vapply(built$x$data, function(tr) {
    length(tr$x %||% tr$y %||% integer(0))
  }, integer(1))

  expect_true(all(trace_lengths > 0),
    info = paste(
      "Empty trace(s) found at index:",
      paste(which(trace_lengths == 0), collapse = ", ")
    )
  )
})

test_that("PercentBarPlot ignores non-numeric experiment_info columns", {
  fx <- make_test_app_with_groups() # has a factor 'cohort' column
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(relativeTo = "none")
      expect_no_warning(output$PercentBar)
    })
  )
})

# ---- Gap 7: flowSOMPieUI / flowSOMPie real rendering -----------------------

test_that("flowSOMPieUI builds a plotOutput once nodes are selected", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(flowSOMPieCols = 4, flowSOMPieSize = 100, flowSOMPieMax = 10)
      ui <- output$flowSOMPieUI
      expect_true(!is.null(ui))
    })
  )
})

test_that("flowSOMPie renders once flowSOMPiePlot() returns a real plot", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(flowSOMPieCols = 4, flowSOMPieSize = 100, flowSOMPieMax = 10)
      out <- output$flowSOMPie
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 8: staticSomGrid with the app's real dList pairs ------------------

test_that("staticSomGrid builds panels/rows for the app's dList", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      session$setInputs(staticSomCols = 2, staticSomHeight = 200)
      ui <- output$staticSomGrid
      expect_true(!is.null(ui))
      # sanity: with 6 default dList pairs and 2 cols, expect >= 1 row structure
    })
  )
})

# ---- Gap 9: downloadPlots filename() and content() -------------------------
#
# NOTE: with the default nPlots = 6, the content() handler's
# `for (plotIdx in seq_len(nPlots)) { dims <- dimSel[[plotIdx]]$dims }` will
# throw "subscript out of bounds" for plotIdx >= 2, because dimSelection()
# is only ever populated with ONE entry (see note above). Using nPlots = 1
# here to get a clean, passing coverage test; consider a *separate* test
# that documents/guards the nPlots > 1 crash as its own bug report:
#
#   test_that("downloadPlots errors with nPlots > 1 (KNOWN ISSUE)", {
#     fx <- make_test_app()  # default nPlots = 6
#     testServer(fx$app, {
#       select_dims(session, "marker1", "marker2")
#       select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
#       session$setInputs(clusterNameSelect = "x")
#       expect_error(session$getDownload("downloadPlots"), "subscript out of bounds")
#     })
#   })

test_that(".writeClusterSelectorPdf writes a non-empty PDF", {
  fx <- make_test_app()
  tmp <- tempfile(fileext = ".pdf")

  # Build a minimal getBasePlotFn that returns a simple ggplot for any dim pair.
  getBasePlotFn <- function(xVar, yVar, colorVar, sizeVar) {
    ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y))
  }

  .writeClusterSelectorPdf(
    file = tmp,
    rs = 1:3,
    sce = fx$sce,
    sceRN = rownames(fx$sce),
    sceCN = colnames(fx$sce),
    metaD = S4Vectors::metadata(fx$sce),
    somCodesName = "SOM_codes",
    dendPlotObj = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:20, 5)))),
    dimPairs = list(c("marker1", "marker2")),
    getBasePlotFn = getBasePlotFn,
    dlColorVar = "n",
    dlSizeVar = "max",
    pctl = 0.99,
    somRasterXy = data.frame(x = 1, y = 1),
    baseRasterGgplot = ggplot2::ggplot(),
    plotRegistry = list()
  )
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})
