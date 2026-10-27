# Tests for the ridge-plot -> main SOM 2D plot link.
#
# Clicking an outlier dot on a ridge row must (a) carry the row's marker in
# the dot's customdata, (b) re-axis the main 2D plot to a marker pair that
# contains that marker -- the first preset pair when one exists, otherwise
# marker-on-x with the first other marker on y -- and (c) highlight the
# clicked cell and its SOM node in the 2D plot.

library(testthat)
library(shiny)

# Local helper: the ones in test-outputs_clusterSelector-gaps.R live in that
# file's environment, so this file carries its own minimal versions.
#
# Note on the dropdown updates: shiny's MockShinySession$sendInputMessage is
# a noop, so the updateSelectInput calls of the observer never appear in
# session$input under testServer. The tests below therefore assert the
# observer's observable side effects (activePlot, dimSelection) and simulate
# the browser's application of the dropdown update with session$setInputs.
ridge_click_json <- function(marker, cell) {
    jsonlite::toJSON(list(list(curveNumber = 0, pointNumber = 0,
        customdata = paste(marker, cell, sep = "\x1f"))), auto_unbox = TRUE)
}

# ---- .ridgeAxisForMarker: the axis decision ---------------------------------

test_that(".ridgeAxisForMarker keeps the axes when the marker is already shown", {
    markers <- c("marker1", "marker2", "marker3")
    dList <- list(c("marker1", "marker2"))
    expect_null(.ridgeAxisForMarker("marker1", c("marker1", "marker2"),
        dList, markers))
    expect_null(.ridgeAxisForMarker("marker2", c("marker1", "marker2"),
        dList, markers))
})

test_that(".ridgeAxisForMarker picks the first preset pair containing the marker", {
    markers <- paste0("marker", seq_len(6))
    dList <- list(c("marker5", "marker6"), c("marker3", "marker4"))
    move <- .ridgeAxisForMarker("marker3", c("marker1", "marker2"), dList,
        markers)
    expect_equal(move$pair, c("marker3", "marker4"))
    expect_equal(move$idx, 2L)
    # The pair keeps its own orientation: the marker may land on y.
    move2 <- .ridgeAxisForMarker("marker4", c("marker1", "marker2"), dList,
        markers)
    expect_equal(move2$pair, c("marker3", "marker4"))
    expect_equal(move2$idx, 2L)
})

test_that(".ridgeAxisForMarker skips degenerate preset pairs", {
    markers <- c("marker1", "marker3", "marker4", "marker5")
    dList <- list("marker5", c("marker3", "marker3"), c("marker4", "marker3"))
    move <- .ridgeAxisForMarker("marker3", c("marker1", "marker2"), dList,
        markers)
    expect_equal(move$pair, c("marker4", "marker3"))
    expect_equal(move$idx, 3L)
})

test_that(".ridgeAxisForMarker falls back to first-marker-on-y", {
    # No preset pair contains the marker: x takes the marker, y takes the
    # first entry of the marker list.
    expect_equal(
        .ridgeAxisForMarker("marker9", c("marker5", "marker6"),
            list(c("marker1", "marker2")), c("marker1", "marker9"))$pair,
        c("marker9", "marker1"))
    # The marker itself is first in the list: y skips it.
    expect_equal(
        .ridgeAxisForMarker("marker1", c("marker2", "marker3"),
            list(c("marker4", "marker5")), c("marker1", "marker2"))$pair,
        c("marker1", "marker2"))
    # No preset pairs at all.
    expect_equal(
        .ridgeAxisForMarker("marker3", c("marker1", "marker2"), list(),
            c("marker1", "marker3"))$pair,
        c("marker3", "marker1"))
    # Degenerate: the marker is the only one in the list -- no change.
    expect_null(.ridgeAxisForMarker("marker1", c("marker2", "marker3"),
        list(), "marker1"))
    # Unknown marker -- no change.
    expect_null(.ridgeAxisForMarker("nope", c("marker1", "marker2"),
        list(c("marker1", "marker2")), c("marker1", "marker2")))
})

# ---- .parseRidgeClick: the click payload ------------------------------------

test_that(".parseRidgeClick reads marker and cell from the dot customdata", {
    d <- data.frame(customdata = I(list("marker3\x1f7")))
    click <- .parseRidgeClick(d)
    expect_equal(click$marker, "marker3")
    expect_equal(click$cell, 7L)
    # plotly may hand back a plain character vector instead of a list column.
    expect_equal(.parseRidgeClick(data.frame(customdata = "marker3\x1f7")),
        click)
    # A bare numeric cell (the "clicked cell" overlay trace) has no marker.
    legacy <- .parseRidgeClick(data.frame(customdata = 12))
    expect_equal(legacy$cell, 12L)
    expect_null(legacy$marker)
    # Malformed or empty payloads parse to NULL.
    expect_null(.parseRidgeClick(data.frame(customdata = "marker3\x1fnope")))
    expect_null(.parseRidgeClick(data.frame(customdata = "nope")))
    expect_null(.parseRidgeClick(NULL))
    expect_null(.parseRidgeClick(data.frame()))
})

test_that(".parseRidgeClick tolerates non-data-frame payloads", {
    # A bare numeric payload (plotly may hand back an atomic value) parses
    # as a legacy cell click.
    expect_equal(.parseRidgeClick(5), list(marker = NULL, cell = 5L))
    # A bare non-numeric payload and NA parse to NULL.
    expect_null(.parseRidgeClick("nope"))
    expect_null(.parseRidgeClick(NA_real_))
    # An empty event ([]) parses to NULL.
    expect_null(.parseRidgeClick(list()))
    # Point objects in plain list form (parse_json without
    # simplifyDataFrame) still carry marker and cell.
    expect_equal(.parseRidgeClick(list(customdata = "marker2\x1f3")),
        list(marker = "marker2", cell = 3L))
    expect_equal(.parseRidgeClick(list(list(customdata = "marker2\x1f3"))),
        list(marker = "marker2", cell = 3L))
})

# ---- plotRidgeInteractive: the dots carry their marker ----------------------

test_that("ridge outlier dots encode their marker in customdata", {
    sce <- make_interactive_sce()
    p <- plotRidgeInteractive(sce, c("b", "a"), rs = c(1L, 2L),
        coverage = 0.9)
    b <- plotly::plotly_build(p)
    dots <- b$x$data[vapply(b$x$data, function(t) {
        grepl("^outlier", if (is.null(t$name)) "" else t$name)
    }, logical(1))]
    expect_true(length(dots) > 0L)
    for (tr in dots) {
        expect_true(is.character(tr$customdata))
        parts <- strsplit(tr$customdata, "\x1f", fixed = TRUE)
        expect_true(all(vapply(parts, function(pp) {
            length(pp) == 2L && pp[1] %in% c("a", "b") &&
                !is.na(suppressWarnings(as.integer(pp[2])))
        }, logical(1))))
    }
})

# ---- Highlight: the clicked cell and its node in the 2D plot ----------------

test_that(".ridgeHighlightPoints locates the clicked cell and its node", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    cell <- 3L
    node <- as.character(sce$cluster_id[cell])
    hl <- .ridgeHighlightPoints(sce, codes, cell, c("b", "a"))
    expect_equal(hl$x, unname(SummarizedExperiment::assay(sce)["b", cell]))
    expect_equal(hl$y, unname(SummarizedExperiment::assay(sce)["a", cell]))
    expect_equal(hl$node, node)
    expect_equal(hl$nodeX, codes[node, "b"])
    expect_equal(hl$nodeY, codes[node, "a"])
})

test_that(".ridgeHighlightPoints falls back to per-node means", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    # Marker "a" is not a SOM column: the node coordinate must be the
    # per-node mean of the cells assigned to the node, the way
    # plotSOMScatter computes it.
    codes2 <- codes[, colnames(codes) != "a", drop = FALSE]
    cell <- 5L
    node <- as.character(sce$cluster_id[cell])
    hl <- .ridgeHighlightPoints(sce, codes2, cell, c("b", "a"))
    expect_equal(hl$nodeY,
        mean(SummarizedExperiment::assay(sce)["a", sce$cluster_id == node]))
    expect_equal(hl$nodeX, codes2[node, "b"])
    expect_equal(hl$node, node)
})

test_that(".ridgeHighlightPoints is safe at the edges", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    expect_null(.ridgeHighlightPoints(sce, codes, NULL, c("b", "a")))
    # Unknown markers: nothing to draw.
    expect_null(.ridgeHighlightPoints(sce, codes, 2L, c("nope", "a")))
    # A cell whose node cannot be resolved still gets cell coordinates.
    sce$cluster_id[2] <- NA
    hl <- .ridgeHighlightPoints(sce, codes, 2L, c("b", "a"))
    expect_equal(hl$x, unname(SummarizedExperiment::assay(sce)["b", 2L]))
    expect_true(is.na(hl$nodeX))
    expect_true(is.na(hl$nodeY))
    expect_true(is.na(hl$node))
})

test_that("highlight_df drops markers missing from the SOM codes", {
    # The axis dropdowns offer all markers, so a pair may contain a marker
    # without a SOM code column. The red node layer must degrade to empty
    # instead of erroring (regression guard for the ridge-click fallback).
    sce <- make_interactive_sce()
    metaD <- S4Vectors::metadata(sce)
    expect_s3_class(
        highlight_df("nope", "b", c(1L, 2L), "SOM_codes", metaD = metaD),
        "data.frame")
    expect_equal(nrow(highlight_df("nope", "b", c(1L, 2L), "SOM_codes",
        metaD = metaD)), 0L)
})

test_that("somPlot draws the clicked cell and its node in the plain view", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    hl <- .ridgeHighlightPoints(sce, codes, 3L, c("b", "a"))
    pp1 <- plotSOMScatter(sce, chs = c("b", "a"))
    p <- somPlot(pp1, 1L, c(1L, 2L), character(0), FALSE,
        dimSelection = list(list(dims = c("b", "a"))), sce = sce,
        metaD = S4Vectors::metadata(sce), cellHighlight = hl)
    b <- plotly::plotly_build(p)
    # The dashed connector between the node and the cell.
    seg <- Filter(function(t) identical(t$mode, "lines"), b$x$data)
    expect_true(any(vapply(seg, function(t) {
        identical(sort(t$x), sort(c(hl$nodeX, hl$x))) &&
            identical(sort(t$y), sort(c(hl$nodeY, hl$y)))
    }, logical(1))))
    # The cell marker is a diamond at the cell's coordinates.
    dia <- Filter(function(t) identical(t$marker$symbol, "diamond"),
        b$x$data)
    expect_true(any(vapply(dia, function(t) {
        length(t$x) == 1L && isTRUE(all.equal(as.numeric(t$x), hl$x)) &&
            isTRUE(all.equal(as.numeric(t$y), hl$y))
    }, logical(1))))
    # The node is ringed (an open circle around the node's coordinates).
    ring <- Filter(function(t) identical(t$marker$symbol, "circle-open"),
        b$x$data)
    expect_true(any(vapply(ring, function(t) {
        length(t$x) == 1L && isTRUE(all.equal(as.numeric(t$x), hl$nodeX)) &&
            isTRUE(all.equal(as.numeric(t$y), hl$nodeY))
    }, logical(1))))
})

test_that("somPlot without a clicked cell draws no highlight layers", {
    sce <- make_interactive_sce()
    pp1 <- plotSOMScatter(sce, chs = c("b", "a"))
    p <- somPlot(pp1, 1L, c(1L, 2L), character(0), FALSE,
        dimSelection = list(list(dims = c("b", "a"))), sce = sce,
        metaD = S4Vectors::metadata(sce))
    b <- plotly::plotly_build(p)
    expect_false(any(vapply(b$x$data, function(t) {
        identical(t$marker$symbol, "diamond")
    }, logical(1))))
})

test_that("a ridge click highlights the cell and node in the main 2D plot", {
    test_app <- make_test_app()
    quiet_plotly_test(
        testServer(app = test_app$app, expr = {
            has_diamond <- function(q) {
                if (is.character(q)) {
                    q <- jsonlite::fromJSON(q, simplifyVector = FALSE)
                }
                any(vapply(q$x$data, function(t) {
                    identical(t$marker$symbol, "diamond")
                }, logical(1)))
            }
            session$setInputs(clusterNumbers = "1,2")
            session$elapse(1001) # inputClusterNumber debounce
            session$elapse(1501) # rsUsed_d debounce
            # The plain (non-projection) view: the checkbox has no UI default
            # under testServer.
            session$setInputs(showGroups = FALSE, nnAgreement = FALSE)
            session$setInputs(currentDimX = "marker1",
                currentDimY = "marker2")
            session$elapse(300)
            expect_false(has_diamond(output$somDataMain))
            session$setInputs(`plotly_click-ridgePlot` =
                ridge_click_json("marker1", 2L))
            session$flushReact()
            after <- output$somDataMain
            expect_true(has_diamond(after))
            # The diamond sits at the cell's values on the active axes.
            if (is.character(after)) {
                after <- jsonlite::fromJSON(after, simplifyVector = FALSE)
            }
            cx <- unname(SummarizedExperiment::assay(
                test_app$sce)["marker1", 2L])
            cy <- unname(SummarizedExperiment::assay(
                test_app$sce)["marker2", 2L])
            dia <- Filter(function(t) identical(t$marker$symbol, "diamond"),
                after$x$data)
            expect_true(any(vapply(dia, function(t) {
                isTRUE(all.equal(unlist(t$x), cx)) &&
                    isTRUE(all.equal(unlist(t$y), cy))
            }, logical(1))))
        })
    )
})

# ---- Layout: the ridge box sits above the SOM 2D plots box ------------------

test_that("the ridge box is rendered above the SOM 2D plots box", {
    app <- make_test_app()
    e <- environment(app$app$serverFuncSource)
    html <- as.character(e$ui)
    ridgePos <- regexpr('id="boxInfo_ridgeBox"', html, fixed = TRUE)
    somPos <- regexpr('id="boxInfo_som2dPlots"', html, fixed = TRUE)
    expect_true(ridgePos > 0L && somPos > 0L)
    expect_lt(ridgePos, somPos)
})

# ---- The ridge-click axis observer (full app) -------------------------------

test_that("a ridge click re-axes the 2D plot to the first matching preset pair", {
    test_app <- make_test_app()
    quiet_plotly_test(
        testServer(app = test_app$app, expr = {
            session$setInputs(currentDimX = "marker1", currentDimY = "marker2")
            session$elapse(300) # activeDims debounce
            session$setInputs(`plotly_click-ridgePlot` =
                ridge_click_json("marker3", 2L))
            session$flushReact()
            # The observer activates the matched preset pair.
            expect_equal(activePlot(), 2L)
            # Simulate the browser applying the observer's dropdown update.
            session$setInputs(currentDimX = "marker3", currentDimY = "marker4")
            session$elapse(300)
            expect_equal(dimSelection()[[1L]]$dims, c("marker3", "marker4"))
        })
    )
})

test_that("a ridge click leaves the preset alone when the marker is already shown", {
    test_app <- make_test_app()
    quiet_plotly_test(
        testServer(app = test_app$app, expr = {
            session$setInputs(currentDimX = "marker1", currentDimY = "marker3")
            session$elapse(300)
            session$setInputs(`plotly_click-ridgePlot` =
                ridge_click_json("marker3", 2L))
            session$flushReact()
            expect_equal(activePlot(), 1L)
            expect_equal(dimSelection()[[1L]]$dims, c("marker1", "marker3"))
        })
    )
})

test_that("a ridge click without a preset pair activates no preset", {
    test_app <- make_test_app()
    quiet_plotly_test(
        testServer(app = test_app$app, expr = {
            dListRV(list(c("marker1", "marker2")))
            session$setInputs(currentDimX = "marker5", currentDimY = "marker6")
            session$elapse(300)
            session$setInputs(`plotly_click-ridgePlot` =
                ridge_click_json("marker9", 2L))
            session$flushReact()
            # No preset pair contains the marker, so no preset is activated;
            # the fallback pair (marker9, marker1) goes to the dropdowns only
            # (the decision itself is unit-tested above).
            expect_equal(activePlot(), 1L)
        })
    )
})

test_that("clicks without a marker payload are ignored without error", {
    test_app <- make_test_app()
    quiet_plotly_test(
        testServer(app = test_app$app, expr = {
            session$setInputs(currentDimX = "marker1", currentDimY = "marker2")
            session$elapse(300)
            # A bare cell index (the "clicked cell" overlay trace) and an
            # empty event must both be no-ops -- in particular they must not
            # raise an error in the observer (shiny would print it).
            for (payload in list("5", "[]")) {
                session$setInputs(`plotly_click-ridgePlot` = payload)
                session$flushReact()
                expect_equal(activePlot(), 1L)
                expect_equal(dimSelection()[[1L]]$dims,
                    c("marker1", "marker2"))
            }
        })
    )
})
