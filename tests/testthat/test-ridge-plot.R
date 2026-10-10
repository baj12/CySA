# tests/testthat/test-ridge-plot.R
#
# Tests for the "marker ridge" box: per-marker ridgeline densities with
# connected trace lines for cells of the current selection that are extreme
# (|value - mean| > threshold * sd) on at least one marker.
#
# Layout: one row per marker (shared x-axis), all-cell density as a filled
# pastel area, selected-cluster density as a solid line, dashed threshold
# ticks at mean +/- threshold * sd, and capped trace lines coloured by the
# marker on which the cell is most extreme.

library(testthat)
library(shiny)

# ---- Local fixture: tiny hand-checkable SCE --------------------------------

# 8 cells, 2 markers, 2 clusters:
#   cluster "1" = cells 1-2, cluster "2" = cells 3-8.
# Values chosen so the extreme-cell rule can be recomputed with base R.
make_ridge_sce <- function() {
    markers <- c("m1", "m2")
    counts <- matrix(
        data = c(
            -5, 0,
            -4, 1,
            0, -1,
            1, 2,
            2, -2,
            3, 3,
            4, -3,
            5, 4
        ),
        nrow = 2, ncol = 8, byrow = TRUE,
        dimnames = list(markers, paste0("cell", seq_len(8)))
    )
    cd <- S4Vectors::DataFrame(
        sample_id = rep(c("sample1", "sample2"), length.out = 8),
        cluster_id = factor(rep(c("1", "1", rep("2", 6))),
            levels = c("1", "2")
        )
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(exprs = counts),
        colData = cd
    )
    som_codes <- matrix(
        data = seq_len(2 * 2) / 4,
        nrow = 2, ncol = 2,
        dimnames = list(c("1", "2"), markers)
    )
    S4Vectors::metadata(sce)$SOM_codes <- som_codes
    S4Vectors::metadata(sce)$SOM_stats <- data.frame(
        id = seq_len(2), n = c(2L, 6L)
    )
    S4Vectors::metadata(sce)$map <- list(colsUsed = markers)
    sce
}

# Base-R reference: which selected-cluster cells are extreme on >= 1 marker?
reference_extreme_cells <- function(sce, markers, rs, threshold) {
    mat <- SummarizedExperiment::assays(sce)[[1]]
    cid <- as.integer(as.character(
        SummarizedExperiment::colData(sce)$cluster_id
    ))
    sel <- which(cid %in% rs)
    if (length(sel) == 0) {
        return(data.frame(
            cellId = integer(0),
            marker = character(0),
            z = numeric(0)
        ))
    }
    stats <- vapply(markers, function(m) {
        v <- mat[m, ]
        c(mean = mean(v), sd = stats::sd(v))
    }, numeric(2))
    extreme <- list()
    for (cell in sel) {
        zs <- vapply(markers, function(m) {
            (mat[m, cell] - stats["mean", m]) / stats["sd", m]
        }, numeric(1))
        hit <- names(zs)[abs(zs) > threshold]
        if (length(hit) > 0) {
            best <- hit[which.max(abs(zs[hit]))]
            extreme[[as.character(cell)]] <- data.frame(
                cellId = cell, marker = best, z = zs[best]
            )
        }
    }
    if (length(extreme) == 0) {
        return(data.frame(
            cellId = integer(0),
            marker = character(0),
            z = numeric(0)
        ))
    }
    do.call(rbind, extreme)
}


test_that(".ridgeAllDensities returns one shared grid per marker with stats", {
    sce <- make_ridge_sce()
    res <- .ridgeAllDensities(sce, "SOM_codes", c("m1", "m2"), nGrid = 64)

    expect_s3_class(res$allDensity, "data.frame")
    expect_true(all(c("marker", "x", "y") %in% colnames(res$allDensity)))
    expect_equal(nrow(res$allDensity[res$allDensity$marker == "m1", ]), 64)
    expect_equal(nrow(res$allDensity[res$allDensity$marker == "m2", ]), 64)
    # Shared x-range: both markers evaluated on the same grid.
    expect_equal(
        sort(unique(res$allDensity$x[res$allDensity$marker == "m1"])),
        sort(unique(res$allDensity$x[res$allDensity$marker == "m2"]))
    )
    expect_true("m1" %in% names(res$markerStats))
    expect_true(all(c("mean", "sd") %in% names(res$markerStats$m1)))
    # Stats match base R.
    mat <- SummarizedExperiment::assays(sce)[[1]]
    expect_equal(unname(res$markerStats$m1["mean"]), mean(mat["m1", ]))
    expect_equal(unname(res$markerStats$m1["sd"]), stats::sd(mat["m1", ]))
})

test_that(".ridgeTraceData extreme-cell rule matches base-R reference", {
    sce <- make_ridge_sce()
    markers <- c("m1", "m2")
    res <- .ridgeTraceData(sce, "SOM_codes", markers,
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    ref <- reference_extreme_cells(sce, markers, rs = c(2), threshold = 0.99)

    got <- unique(res$traces$cellId[res$traces$isExtreme])
    expect_setequal(got, ref$cellId)
})

test_that(".ridgeTraceData keeps only the most extreme cells under the cap", {
    sce <- make_ridge_sce()
    markers <- c("m1", "m2")
    res <- .ridgeTraceData(sce, "SOM_codes", markers,
        rs = c(1, 2), threshold = 0.99, maxLines = 2
    )
    ref <- reference_extreme_cells(sce, markers, rs = c(1, 2), threshold = 0.99)
    ref <- ref[order(-abs(ref$z)), ]

    got <- unique(res$traces$cellId[res$traces$isExtreme])
    expect_length(got, 2)
    expect_setequal(got, ref$cellId[seq_len(2)])
})

test_that(".ridgeTraceData marks extreme rows and colours by extreme marker", {
    sce <- make_ridge_sce()
    markers <- c("m1", "m2")
    res <- .ridgeTraceData(sce, "SOM_codes", markers,
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    ref <- reference_extreme_cells(sce, markers, rs = c(2), threshold = 0.99)

    tr <- res$traces[res$traces$cellId %in% ref$cellId, ]
    expect_true(all(tr$isExtreme[tr$marker == ref$marker[1] &
        tr$cellId == ref$cellId[1]]))
    # Trace rows cover every selected marker, not only the extreme one.
    expect_setequal(
        unique(res$traces$marker[res$traces$cellId == ref$cellId[1]]),
        markers
    )
    # Extreme marker recorded for colouring.
    expect_true(all(
        res$traceMeta$extremeMarker[res$traceMeta$cellId == ref$cellId[1]]
        == ref$marker[1]
    ))
})

test_that(".ridgeTraceData returns empty plot data for empty selection", {
    sce <- make_ridge_sce()
    res <- .ridgeTraceData(sce, "SOM_codes", c("m1", "m2"),
        rs = integer(0), threshold = 0.99, maxLines = 50
    )
    expect_null(res$traces)
    expect_null(res$traceMeta)
    # The all-cell densities do not depend on the selection: the same
    # densities come back regardless of rs.
    d_with <- .ridgeAllDensities(sce, "SOM_codes", c("m1", "m2"), nGrid = 8)
    d_without <- .ridgeAllDensities(sce, "SOM_codes", c("m1", "m2"), nGrid = 8)
    expect_equal(d_with$allDensity, d_without$allDensity)
})

test_that("trace data tolerates an assay without cell colnames", {
    # Real CATALYST SCEs can lose assay colnames (e.g. when round-tripped
    # through saveRDS of an empty-assay subsample); the trace extraction
    # must index positionally instead of by cell name.
    sce <- make_ridge_sce()
    mat <- SummarizedExperiment::assays(sce)[[1]]
    naked <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = matrix(mat,
            nrow = nrow(mat),
            dimnames = list(rownames(mat), NULL)
        )),
        colData = SummarizedExperiment::colData(sce)
    )
    S4Vectors::metadata(naked) <- S4Vectors::metadata(sce)
    expect_no_error(res <- .ridgeTraceData(naked, "SOM_codes",
        c("m1", "m2"), rs = c(2), threshold = 0.99, maxLines = 50
    ))
    ref <- reference_extreme_cells(sce, c("m1", "m2"), rs = c(2),
        threshold = 0.99
    )
    got <- unique(res$traces$cellId[res$traces$isExtreme])
    expect_setequal(got, ref$cellId)
    # The traced values match the named-colnames run.
    withNames <- .ridgeTraceData(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    expect_equal(res$traces$x[order(res$traces$cellId, res$traces$marker)],
        withNames$traces$x[order(withNames$traces$cellId,
            withNames$traces$marker)]
    )
})

test_that("trace lines are connected across marker rows in one panel", {
    sce <- make_ridge_sce()
    markers <- c("m1", "m2")
    p <- plotRidgeTraceFunc(sce, "SOM_codes", markers,
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    expect_length(tracedIdx, 1)
    b <- ggplot2::ggplot_build(p)
    d <- b$data[[tracedIdx]]
    d <- d[d$datatype == "path", ]
    # The traces are drawn in a single coordinate system (one panel), so a
    # cell's path actually connects its values across the marker rows.
    expect_length(unique(d$PANEL), 1)
    # Each traced cell covers every marker row.
    expect_true(all(table(d$group) == length(markers)))
})

test_that("plotRidgeTraceFunc returns a ggplot with ridge layers", {
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    expect_s3_class(p, "ggplot")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    # The ridges, traces and dots live in the single traced layer.
    expect_true("GeomRidgelineTraced" %in% geoms)
})

test_that("ridges, traces and dots render in one GeomRidgelineTraced layer", {
    # The ggridges-based rework renders the ridges (ridgeline datatype),
    # the connected trace polylines (path datatype) and the outlier dots
    # (point datatype) in a single GeomRidgelineTraced layer, so the
    # trace z-order (above ridges, below dots) is guaranteed by the geom.
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    traced <- which(geoms == "GeomRidgelineTraced")
    expect_length(traced, 1)
    ld <- ggplot2::ggplot_build(p)$data[[traced]]
    # All three datatypes present in the one layer.
    expect_setequal(unique(ld$datatype),
        c("ridgeline", "path", "point"))
    # The traces connect across rows: the path group of a traced cell has
    # one vertex per marker row.
    pathGroups <- unique(ld$group[ld$datatype == "path"])
    expect_true(length(pathGroups) > 0)
    expect_true(all(table(ld$group[ld$datatype == "path"]) ==
        length(unique(ld$marker[ld$datatype == "ridgeline"]))))
    # Dots are drawn for extreme vertices only, grouped with the cell's
    # path.
    dotGroups <- unique(ld$group[ld$datatype == "point"])
    expect_true(all(dotGroups %in% pathGroups))
    # The trace group's ymin is below every ridge base, so the draw_panel
    # ordering puts it on top of all ridges.
    expect_lt(max(ld$ymin[ld$datatype %in% c("path", "point")]),
        min(ld$ymin[ld$datatype == "ridgeline"]))
    # Dots sit on their cluster's ridge: at least one above its row base
    # by more than rounding.
    expect_true(any(ld$ymax[ld$datatype == "point"] >
        ld$ymin[ld$datatype == "point"] + 1e-4))
})

test_that("showTraces=FALSE drops path rows but keeps point rows", {
    # The single-layer rework: showTraces=FALSE drops only the "path"
    # datatype rows; the outlier dots ("point") remain, and the ridges
    # stay.
    sce <- make_ridge_sce()
    pOff <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50, showTraces = FALSE
    )
    geoms <- vapply(pOff$layers, function(l) class(l$geom)[1], character(1))
    traced <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(pOff)$data[[traced]]
    expect_false("path" %in% ld$datatype)
    expect_true("point" %in% ld$datatype)
    expect_true("ridgeline" %in% ld$datatype)
})

test_that("trace vertices sit on the density curve, not flattened at the base", {
    # Bug: trace y was looked up by exact grid match, but a cell's value
    # almost never falls exactly on the density grid, so every vertex
    # collapsed to the row base and the plot looked like a plain ridge
    # plot with no visible outliers. Vertices must be interpolated at the
    # cell's value and sit above the row base.
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    paths <- ld[ld$datatype == "path", ]
    # Every vertex must sit strictly above the bottom of the panel: the
    # flattened bug put vertices at y = 0 within-row (i.e. only the very
    # bottom row had correct y, all others sat at their base exactly).
    # With correct interpolation, each vertex sits at its cluster ridge
    # (row base + cluster offset + density >= row base).
    ridges <- ld[ld$datatype == "ridgeline", ]
    for (r in seq_len(nrow(paths))) {
        m <- as.character(paths$marker[r])
        rowBase <- min(ridges$ymin[as.character(ridges$marker) == m])
        expect_gte(paths$y[r], rowBase - 1e-6)
    }
    # And at least one vertex must sit strictly above its row base
    # (the density at the cell value is positive).
    expect_true(any(vapply(seq_len(nrow(paths)), function(r) {
        m <- as.character(paths$marker[r])
        rowBase <- min(ridges$ymin[as.character(ridges$marker) == m])
        paths$y[r] > rowBase + 1e-4
    }, logical(1))))
})

test_that("plotRidgeTraceFunc handles empty selection and bad markers", {
    sce <- make_ridge_sce()

    p1 <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = integer(0), threshold = 0.99, maxLines = 50
    )
    expect_s3_class(p1, "ggplot")
    expect_true(grepl("no nodes selected", p1$labels$title, ignore.case = TRUE))

    p2 <- plotRidgeTraceFunc(sce, "SOM_codes", c("not_a_marker"),
        rs = c(1), threshold = 0.99, maxLines = 50
    )
    expect_s3_class(p2, "ggplot")
    expect_true(grepl("no valid markers", p2$labels$title, ignore.case = TRUE))
})

test_that("ridge box exposes inputs and output in the UI", {
    test_app <- make_test_app()
    sce <- test_app$sce
    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed
    box <- .buildRidgeBox(som_codes, colsUsed)
    html <- as.character(box)
    expect_true(grepl("ridgeMarkerSelection", html, fixed = TRUE))
    expect_true(grepl("ridgeThreshold", html, fixed = TRUE))
    expect_true(grepl("ridgeMaxLines", html, fixed = TRUE))
    expect_true(grepl('id="ridgePlot"', html, fixed = TRUE))
})

test_that("server ridge output renders and reacts to inputs", {
    test_app <- make_test_app()
    shiny::testServer(
        app = test_app$app,
        expr = {
            suppressWarnings(session$setInputs(
                ridgeMarkerSelection = c("marker1", "marker2")
            ))
            suppressWarnings(session$flushReact())
            expect_true(!is.null(output$ridgePlot))

            suppressWarnings(session$setInputs(
                ridgeThreshold = 0.5, ridgeMaxLines = 5
            ))
            suppressWarnings(session$flushReact())
            expect_true(!is.null(output$ridgePlot))
        }
    )
})

test_that("ridge box renders a cluster picker driven by the selection", {
    test_app <- make_test_app()
    shiny::testServer(
        app = test_app$app,
        expr = {
            # More than 5 selected nodes -> the picker must appear.
            suppressWarnings(session$setInputs(
                clusterNumbers = paste(seq_len(7), collapse = ",")
            ))
            session$flushReact()
            session$elapse(1100)
            suppressWarnings(session$flushReact())
            html <- paste(as.character(output$ridgeClusterUI),
                collapse = " ")
            expect_true(grepl("ridgeClusterSelection", html, fixed = TRUE))
            # All 7 selected nodes are choices; the 5 largest (ties: lowest
            # ids) are pre-selected.
            expect_true(grepl("<option value=\"1\" selected>", html,
                fixed = TRUE))
            expect_true(grepl("<option value=\"7\">", html, fixed = TRUE))
            expect_false(grepl("<option value=\"7\" selected>", html,
                fixed = TRUE))
        }
    )
})

test_that("ridge picker is hidden when 5 or fewer clusters are selected", {
    test_app <- make_test_app()
    shiny::testServer(
        app = test_app$app,
        expr = {
            suppressWarnings(session$setInputs(
                clusterNumbers = paste(seq_len(4), collapse = ",")
            ))
            session$flushReact()
            session$elapse(1100)
            suppressWarnings(session$flushReact())
            expect_null(output$ridgeClusterUI)
        }
    )
})

test_that("server ridge honours the cluster picker when it is set", {
    test_app <- make_test_app()
    shiny::testServer(
        app = test_app$app,
        expr = {
            suppressWarnings(session$setInputs(
                clusterNumbers = paste(seq_len(7), collapse = ",")
            ))
            session$flushReact()
            session$elapse(1100)
            suppressWarnings(session$flushReact())
            suppressWarnings(session$setInputs(
                ridgeClusterSelection = "1"
            ))
            suppressWarnings(session$flushReact())
            session$elapse(600) # ridgePlot debounce (500 ms)
            suppressWarnings(session$flushReact())
            ridgeServer <- get("ridgeServer",
                envir = environment(session$testserver_serverfn)
            )
            b <- plotly::plotly_build(ridgeServer$ridgePlot())
            nms <- unique(vapply(b$x$data, function(t) {
                if (is.null(t$name)) "" else t$name
            }, character(1)))
            expect_true("node 1" %in% nms)
            expect_false(any(sprintf("node %d", 2:7) %in% nms))
            expect_true("selection" %in% nms)
        }
    )
})

test_that("rows are separated with headroom and no overlapping bands", {
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    # The traced layer's ridgeline rows carry the row offsets with
    # headroom: the base of each row must sit above the peak of the row
    # below.
    ld <- ggplot2::ggplot_build(p)$data[[
        which(vapply(p$layers, function(l) class(l$geom)[1],
            character(1)) == "GeomRidgelineTraced")]]
    ridges <- ld[ld$datatype == "ridgeline", ]
    # The top row's base must sit above the lower row's peak, so no band
    # reaches into the row above it (they are visually separated). rowHeight
    # accounts for the cluster-stack height, so the plain separation holds.
    lowerPeak <- max(ridges$ymax[ridges$marker == "m2"])
    upperBase <- min(ridges$ymin[ridges$marker == "m1"])
    expect_gt(upperBase, lowerPeak)
})

test_that("thresholds are dashed vertical lines spanning each row", {
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomSegment" %in% geoms)
    segIdx <- which(geoms == "GeomSegment")
    seg <- p$layers[[segIdx]]$data
    # Two segments per marker (low and high threshold).
    expect_equal(nrow(seg), 2 * length(unique(seg$marker)))
    # Each segment spans from its row base to the row top (headroom).
    expect_true(all(seg$yend > seg$y))
    # Linetype is dashed and linewidth is set.
    pars <- p$layers[[segIdx]]$aes_params
    expect_identical(pars$linetype, "dashed")
    expect_true(is.numeric(pars$linewidth))
})

test_that("traced cells are drawn as dots at their values", {
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    dots <- ld[ld$datatype == "point", ]
    expect_true(nrow(dots) > 0)
    # One dot per traced cell per marker row (markers that are extreme for
    # any traced cell).
    expect_equal(nrow(dots), length(unique(dots$cellId)) *
        length(unique(dots$marker)))
})

test_that("server ridge reactive is debounced against resize storms", {
    # The heavy ridge computation must be debounced so the stream of
    # jqui_resizable resize events does not re-plot continuously. shiny::
    # debounce wraps the reactive and its result carries the
    # "reactive.event" class; a plain reactive does not.
    test_app <- make_test_app()
    shiny::testServer(
        app = test_app$app,
        expr = {
            suppressWarnings(session$setInputs(
                ridgeMarkerSelection = c("marker1", "marker2")
            ))
            suppressWarnings(session$flushReact())
            expect_true(!is.null(output$ridgePlot))
            ridgeServer <- get("ridgeServer",
                envir = environment(session$testserver_serverfn)
            )
            expect_true("reactive.event" %in% class(ridgeServer$ridgePlot))
        }
    )
})

test_that("dashboard header shows the CySA version", {
    # The header title must contain the package version so users and bug
    # reports always reference the right build. shiny::shinyApp() keeps the
    # ui only inside its closure, so call .buildClusterSelectorApp (the
    # path clusterSelector() takes) with the make_test_app fixture pieces
    # and inspect the HTML.
    set.seed(12)
    sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
    S4Vectors::metadata(sce)$experiment_info$treatment <- factor(
        rep(c("control", "treated"),
            length.out = nrow(S4Vectors::metadata(sce)$experiment_info)),
        levels = c("control", "treated")
    )
    prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    markers <- S4Vectors::metadata(sce)$map$colsUsed
    dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
    dendTable <- data.frame(
        id = seq_len(nrow(som_codes)),
        label = rownames(som_codes),
        stringsAsFactors = FALSE
    )
    cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
    somRasterData <- data.frame(
        x = rep(seq_len(5), length.out = nrow(som_codes)),
        y = rep(seq_len(2), each = nrow(som_codes) / 2),
        id = seq_len(nrow(som_codes))
    )
    for (m in markers) {
        somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
    }
    arr <- array(
        seq_len(10 * 10 * length(markers)),
        dim = c(10, 10, length(markers))
    )
    somRasterObj <- raster::brick(arr)
    names(somRasterObj) <- markers
    csa <- CySA:::.buildClusterSelectorApp(
        sce = prepped$sce,
        sce_subsampled = prepped$sce_subsampled,
        dList = prepped$dList,
        dend = dend,
        dendTable = dendTable,
        clusterPatientTable = cpt,
        somRasterData = somRasterData,
        somRasterObj = somRasterObj,
        env = new.env()
    )
    html <- as.character(csa$ui)
    expect_true(grepl("Cluster Selector", html, fixed = TRUE))
    expect_true(grepl("CySA [0-9]+[.][0-9]+[.][0-9]+", html))
})

test_that("trace data records each traced cell's cluster", {
    sce <- make_ridge_sce()
    markers <- c("m1", "m2")
    res <- .ridgeTraceData(sce, "SOM_codes", markers,
        rs = c(1, 2), threshold = 0.99, maxLines = 50
    )
    # traceMeta carries the cluster of each traced cell.
    expect_true("cluster" %in% colnames(res$traceMeta))
    ref <- reference_extreme_cells(sce, markers, rs = c(1, 2), threshold = 0.99)
    for (k in seq_len(nrow(res$traceMeta))) {
        cid <- res$traceMeta$cellId[k]
        expect_identical(res$traceMeta$cluster[k],
            as.character(SummarizedExperiment::colData(sce)$cluster_id[cid]))
    }
})

test_that("per-cluster lines and combined line are drawn per marker row", {
    sce <- make_ridge_sce()
    # clusters 1 and 2 both selected; each gets its own line, plus a
    # combined (all selected) line, plus the all-cell background area.
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(1, 2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    lines <- ld[ld$datatype == "ridgeline" & !is.na(ld$cluster), ]
    # One data frame whose rows cover the combined ridge line plus the two
    # per-cluster ridge lines, each tagged with its cluster id
    # ("combined" for the pooled line).
    expect_setequal(unique(as.character(lines$cluster)),
        c("combined", "1", "2"))
    # The combined line must be thicker so it reads as the aggregate.
    expect_gt(unique(lines$linewidth[lines$cluster == "combined"])[1],
        max(lines$linewidth[lines$cluster != "combined"]))
})

test_that("per-cluster cap falls back to combined when more than 5", {
    # 6 selected clusters with a cap of 5 -> only the combined line.
    sce <- make_ridge_sce()
    # Build a bigger fixture: 7 clusters over 14 cells.
    counts <- matrix(rnorm(2 * 14), nrow = 2,
        dimnames = list(c("m1", "m2"), paste0("cell", seq_len(14))))
    cd <- S4Vectors::DataFrame(
        sample_id = rep(c("s1", "s2"), length.out = 14),
        cluster_id = factor(rep(seq_len(7), length.out = 14),
            levels = as.character(seq_len(7)))
    )
    sce7 <- SingleCellExperiment::SingleCellExperiment(
        assays = list(exprs = counts), colData = cd)
    S4Vectors::metadata(sce7)$SOM_codes <- matrix(0.5, 7, 2,
        dimnames = list(as.character(seq_len(7)), c("m1", "m2")))
    S4Vectors::metadata(sce7)$map <- list(colsUsed = c("m1", "m2"))

    p <- plotRidgeTraceFunc(sce7, "SOM_codes", c("m1", "m2"),
        rs = seq_len(7), threshold = 0.99, maxLines = 50, maxClusters = 5
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    lines <- ld[ld$datatype == "ridgeline" & !is.na(ld$cluster), ]
    expect_setequal(unique(as.character(lines$cluster)), "combined")
})

test_that("traces are coloured by the traced cell's cluster", {
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(1, 2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    paths <- ld[ld$datatype == "path", ]
    # Path data carries the cluster of each cell, used for colouring.
    expect_true("cluster" %in% colnames(paths))
    metaRef <- .ridgeTraceData(sce, "SOM_codes", c("m1", "m2"),
        rs = c(1, 2), threshold = 0.99, maxLines = 50
    )
    for (k in seq_len(nrow(metaRef$traceMeta))) {
        cid <- metaRef$traceMeta$cellId[k]
        expect_identical(
            unique(as.character(paths$cluster[paths$cellId == cid])),
            metaRef$traceMeta$cluster[k]
        )
    }
})

test_that("cluster line colours are distinct from the background palette", {
    # Bug: colorRampPalette(Set2)(n) interpolates across the whole ramp, so
    # with 2 clusters the second colour is a washed-out grey that is
    # indistinguishable from the default dropped-colour grey and from the
    # marker fill family. Distinct clusters must map to distinct,
    # saturated brewer colours so they read as separate entities.
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(1, 2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    ridges <- ld[ld$datatype == "ridgeline" & !is.na(ld$cluster) &
        ld$cluster != "combined", ]
    vals <- unique(ridges$colour)
    # Two different clusters -> two different, non-grey colours.
    expect_length(vals, 2)
    expect_false(isTRUE(all.equal(vals[1], vals[2])))
    expect_false(identical(unname(vals[2]), "#B3B3B3"))
    # Both must be real hex colours (not NA).
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", vals)))
})

test_that("version tag sits in its own styled line in the header", {
    # Bug: the version was appended inline to "Cluster Selector" inside the
    # 230px-wide AdminLTE logo box with overflow hidden, so long titles
    # clipped it off and users never saw it. It must render on its own line
    # (display:block) so it is always visible.
    set.seed(12)
    sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
    S4Vectors::metadata(sce)$experiment_info$treatment <- factor(
        rep(c("control", "treated"),
            length.out = nrow(S4Vectors::metadata(sce)$experiment_info)),
        levels = c("control", "treated")
    )
    prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    markers <- S4Vectors::metadata(sce)$map$colsUsed
    dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
    dendTable <- data.frame(
        id = seq_len(nrow(som_codes)),
        label = rownames(som_codes),
        stringsAsFactors = FALSE
    )
    cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
    somRasterData <- data.frame(
        x = rep(seq_len(5), length.out = nrow(som_codes)),
        y = rep(seq_len(2), each = nrow(som_codes) / 2),
        id = seq_len(nrow(som_codes))
    )
    for (m in markers) {
        somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
    }
    arr <- array(
        seq_len(10 * 10 * length(markers)),
        dim = c(10, 10, length(markers))
    )
    somRasterObj <- raster::brick(arr)
    names(somRasterObj) <- markers
    csa <- CySA:::.buildClusterSelectorApp(
        sce = prepped$sce,
        sce_subsampled = prepped$sce_subsampled,
        dList = prepped$dList,
        dend = dend,
        dendTable = dendTable,
        clusterPatientTable = cpt,
        somRasterData = somRasterData,
        somRasterObj = somRasterObj,
        env = new.env()
    )
    html <- as.character(csa$ui)
    i <- regexpr("Cluster Selector", html, fixed = TRUE)
    m <- substr(html, i, min(nchar(html), i + 600))
    # The version span must carry display:block (own line) so it cannot be
    # clipped by the fixed-width logo box.
    expect_true(grepl(
        "display\\s*:\\s*block", m
    ))
    # display:block alone is not enough: the logo box is height:50px with
    # line-height:50px inherited by the first text line, which consumes the
    # whole box and pushes the version line below the bottom edge where
    # overflow:hidden clips it vertically. The header must carry a style
    # rule that relaxes the logo line-height so both lines fit in 50px.
    expect_true(grepl(
        "\\.main-header \\.logo\\s*\\{[^}]*line-height\\s*:\\s*(normal|1\\.[0-9]+)",
        html
    ))
})

test_that("per-cluster densities render as stacked joyplot ridges", {
    # Style classique : dans chaque rangée de marqueur, chaque cluster
    # sélectionné est une crête remplie semi-transparente, empilée
    # verticalement (offset par cluster), et non une simple ligne plate.
    sce <- make_ridge_sce()
    p <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(1, 2), threshold = 0.99, maxLines = 50
    )
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    tracedIdx <- which(geoms == "GeomRidgelineTraced")
    ld <- ggplot2::ggplot_build(p)$data[[tracedIdx]]
    areas <- ld[ld$datatype == "ridgeline", ]
    # The combined band and the two per-cluster ridges, each tagged with
    # its cluster id in the cluster column (NA for the background).
    expect_setequal(
        unique(as.character(areas$cluster[!is.na(areas$cluster)])),
        c("combined", "1", "2")
    )
    # Each cluster is offset vertically: same marker, different ridge base
    # between cluster 1 and cluster 2.
    m1 <- areas[areas$marker == "m1" & !is.na(areas$cluster), ]
    base1 <- min(m1$ymin[m1$cluster == "1"])
    base2 <- min(m1$ymin[m1$cluster == "2"])
    expect_false(isTRUE(all.equal(base1, base2)))
})

test_that("plot outputs inside resizable boxes are height-pinned", {
    # Bug: shinyjqui::jqui_resizable wraps each output in a .jqui-wrapper
    # sized with width/height 100% when the output is hidden at wrap time
    # (collapsed boxes, e.g. "marker ridge"). 100% of an auto-height parent
    # is content-driven: the rendered <img> sizes the wrapper, the wrapper
    # sizes the next render, and each render re-reports a new size to the
    # server - an endless render loop on any resize. The header must carry
    # CSS that pins these outputs to a fixed pixel height so the size
    # chain is broken, and the window-resize script must not re-impose
    # 100% heights on image outputs.
    set.seed(12)
    sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
    S4Vectors::metadata(sce)$experiment_info$treatment <- factor(
        rep(c("control", "treated"),
            length.out = nrow(S4Vectors::metadata(sce)$experiment_info)),
        levels = c("control", "treated")
    )
    prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    markers <- S4Vectors::metadata(sce)$map$colsUsed
    dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
    dendTable <- data.frame(
        id = seq_len(nrow(som_codes)),
        label = rownames(som_codes),
        stringsAsFactors = FALSE
    )
    cpt <- table(sample_id = sce$sample_id, cluster_id = sce$cluster_id)
    somRasterData <- data.frame(
        x = rep(seq_len(5), length.out = nrow(som_codes)),
        y = rep(seq_len(2), each = nrow(som_codes) / 2),
        id = seq_len(nrow(som_codes))
    )
    for (m in markers) {
        somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
    }
    arr <- array(
        seq_len(10 * 10 * length(markers)),
        dim = c(10, 10, length(markers))
    )
    somRasterObj <- raster::brick(arr)
    names(somRasterObj) <- markers
    csa <- CySA:::.buildClusterSelectorApp(
        sce = prepped$sce,
        sce_subsampled = prepped$sce_subsampled,
        dList = prepped$dList,
        dend = dend,
        dendTable = dendTable,
        clusterPatientTable = cpt,
        somRasterData = somRasterData,
        somRasterObj = somRasterObj,
        env = new.env()
    )
    html <- as.character(csa$ui)
    # A CSS rule must pin the default plot output height inside boxes.
    expect_true(grepl(
        "\\.tab-content \\.shiny-plot-output\\s*\\{[^}]*height\\s*:\\s*[0-9]+px",
        html
    ))
    # The window-resize script must not re-impose 100% heights on image
    # outputs (the old script did, restarting the content-driven loop).
    expect_false(grepl(
        "ui-resizable[^\\}]*height[^\\}]*100%",
        html
    ))
})

test_that("showTraces=FALSE drops trace lines but keeps outlier dots", {
    # User request: the connecting trace lines made the joyplot unreadable.
    # plotRidgeTraceFunc(showTraces = FALSE) must drop the "path" datatype
    # rows (the polylines across rows) but keep the "point" outlier dots,
    # and default to TRUE so nothing changes for other callers.
    sce <- make_ridge_sce()
    pOff <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50, showTraces = FALSE
    )
    geoms <- vapply(pOff$layers, function(l) class(l$geom)[1], character(1))
    ld <- ggplot2::ggplot_build(pOff)$data[[which(geoms ==
        "GeomRidgelineTraced")]]
    expect_false("path" %in% ld$datatype)
    expect_true("point" %in% ld$datatype)
    # The dots must still sit on the cluster ridge: at least one above
    # its marker row base.
    dots <- ld[ld$datatype == "point", ]
    ridges <- ld[ld$datatype == "ridgeline", ]
    above <- vapply(seq_len(nrow(dots)), function(r) {
        m <- as.character(dots$marker[r])
        rowBase <- min(ridges$ymin[as.character(ridges$marker) == m])
        dots$y[r] > rowBase + 1e-4
    }, logical(1))
    expect_true(any(above))
    # Default keeps the lines.
    pOn <- plotRidgeTraceFunc(sce, "SOM_codes", c("m1", "m2"),
        rs = c(2), threshold = 0.99, maxLines = 50
    )
    ldOn <- ggplot2::ggplot_build(pOn)$data[[which(vapply(pOn$layers,
        function(l) class(l$geom)[1], character(1)) ==
        "GeomRidgelineTraced")]]
    expect_true("path" %in% ldOn$datatype)
})

test_that("trace rows do not stretch the y axis below the ridges", {
    # Regression: the geom gave trace rows a fixed ymin offset of -1 data
    # unit. With expression values in the hundreds the densities are
    # ~1e-3, so the y range was dominated by that offset and every ridge
    # collapsed into a flat line at the top of the panel.
    sce <- CySA_example_sce(n_cells = 500, n_nodes = 25)
    for (st in c(TRUE, FALSE)) {
        p <- plotRidgeTraceFunc(sce, markers = rownames(sce)[seq_len(4)],
            rs = c(1L, 2L, 3L), showTraces = st
        )
        b <- ggplot2::ggplot_build(p)
        ld <- b$data[[1]]
        expect_true(any(ld$datatype == "point"))
        ridges <- ld[ld$datatype == "ridgeline", ]
        ridgeSpan <- diff(range(c(ridges$ymin, ridges$ymax)))
        yRange <- b$layout$panel_params[[1]]$y.range
        # The ridges must fill most of the panel height.
        expect_gt(ridgeSpan / diff(yRange), 0.8)
    }
})

test_that("trace values belong to the traced cell when few cells are extreme", {
    # Regression: trace values were looked up by the cell's position among
    # the extreme cells, but the value matrix is ordered by all selected
    # cells, so traced cells got other cells' values whenever not every
    # selected cell was extreme.
    set.seed(5)
    mat <- rbind(m1 = rnorm(40), m2 = rnorm(40))
    mat[, 5] <- c(9, 0.1)   # only a few selected cells are extreme
    mat[, 12] <- c(-8, 0.2)
    colnames(mat) <- paste0("c", seq_len(40))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(exprs = mat),
        colData = S4Vectors::DataFrame(
            sample_id = "s1",
            cluster_id = factor(rep(c("1", "2"), each = 20))
        )
    )
    tr <- .ridgeTraceData(sce, "SOM_codes", c("m1", "m2"), rs = 1L,
        threshold = 3, maxLines = 50L)$traces
    expect_setequal(unique(tr$cellId), c(5L, 12L))
    expect_equal(tr$x,
        mat[cbind(match(as.character(tr$marker), rownames(mat)), tr$cellId)])
})
