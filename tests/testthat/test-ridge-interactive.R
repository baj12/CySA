# Tests for the interactive marker ridge box (plot-ridge-interactive.R,
# server-ridge.R).

# Two markers on different scales, a selection of nodes 1 and 2; marker
# "b" is listed first in the assay but separates less, so sorting reorders
# the rows relative to the matrix.
make_interactive_sce <- function() {
    set.seed(3)
    n <- 200
    cl <- rep(c(1L, 2L, 3L, 4L), each = n / 4)
    a <- rnorm(n, mean = ifelse(cl <= 2L, 3, 0), sd = 0.5)
    b <- rnorm(n, mean = ifelse(cl <= 2L, 0.2, 0), sd = 1)
    a[1] <- 8 # extreme selected cell on marker a
    mat <- rbind(b = b, a = a)
    colnames(mat) <- paste0("c", seq_len(n))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(exprs = mat),
        colData = S4Vectors::DataFrame(
            sample_id = rep(c("s1", "s2"), length.out = n),
            cluster_id = factor(cl)
        )
    )
    codes <- t(vapply(seq_len(4), function(k) {
        rowMeans(mat[, cl == k, drop = FALSE])
    }, numeric(2)))
    rownames(codes) <- as.character(seq_len(4))
    S4Vectors::metadata(sce)$SOM_codes <- codes
    S4Vectors::metadata(sce)$map <- list(colsUsed = c("b", "a"))
    sce
}

test_that(".ridgeSelectionStats matches base-R means, sds and quantile band", {
    sce <- make_interactive_sce()
    mat <- SummarizedExperiment::assay(sce)
    sel <- which(sce$cluster_id %in% c("1", "2"))
    st <- .ridgeSelectionStats(mat, sel, coverage = 0.9)
    expect_equal(st$marker, c("b", "a"))
    expect_equal(st$selMean, unname(rowMeans(mat[, sel])))
    expect_equal(st$selSd, unname(apply(mat[, sel], 1, sd)))
    q <- apply(mat[, sel], 1, quantile, probs = c(0.05, 0.95),
        names = FALSE)
    expect_equal(st$lower, unname(q[1, ]))
    expect_equal(st$upper, unname(q[2, ]))
    expect_equal(st$separation,
        (st$selMean - st$allMean) / st$allSd)
})

test_that(".ridgeOutliers flags the selected cells outside the band", {
    sce <- make_interactive_sce()
    mat <- SummarizedExperiment::assay(sce)
    sel <- which(sce$cluster_id %in% c("1", "2"))
    st <- .ridgeSelectionStats(mat, sel, coverage = 0.9)
    o <- .ridgeOutliers(mat, sel, st, maxCells = 1000L)
    expect_true(all(o$cell %in% sel))
    # Reference: every selected (cell, marker) outside the quantile band.
    ref <- which(mat[, sel] < st$lower | mat[, sel] > st$upper,
        arr.ind = TRUE)
    expect_equal(nrow(o), nrow(ref))
    # About 1 - coverage of the selected cells fall outside per marker.
    nSel <- length(sel)
    for (m in rownames(mat)) {
        expect_lte(abs(sum(o$marker == m) / nSel - 0.1), 0.02)
    }
    # Values are the real assay values of that cell and marker.
    expect_equal(o$value, mat[cbind(match(o$marker, rownames(mat)), o$cell)])
    expect_true(1L %in% o$cell[o$marker == "a"])
})

test_that(".ridgeOutliers reads values by marker name when rows are sorted", {
    # Regression: stats sorted by separation but mat in assay order used to
    # attach the wrong marker's values to the dots.
    sce <- make_interactive_sce()
    mat <- SummarizedExperiment::assay(sce)
    sel <- which(sce$cluster_id %in% c("1", "2"))
    st <- .ridgeSelectionStats(mat, sel, coverage = 0.9)
    st <- st[order(-abs(st$separation)), ]
    expect_equal(st$marker, c("a", "b"))
    o <- .ridgeOutliers(mat, sel, st, maxCells = 1000L)
    expect_equal(o$value, mat[cbind(match(o$marker, rownames(mat)), o$cell)])
    zRef <- (o$value - st$selMean[match(o$marker, st$marker)]) /
        st$selSd[match(o$marker, st$marker)]
    expect_equal(o$z, zRef)
})

test_that(".ridgeOutliers caps to the most extreme cells", {
    sce <- make_interactive_sce()
    mat <- SummarizedExperiment::assay(sce)
    sel <- which(sce$cluster_id %in% c("1", "2"))
    st <- .ridgeSelectionStats(mat, sel, coverage = 0.5)
    o <- .ridgeOutliers(mat, sel, st, maxCells = 1L)
    expect_equal(unique(o$cell), 1L)
    expect_equal(nrow(.ridgeOutliers(mat, sel[1], st)), 0L)
})

test_that(".cellMembership decomposes the distance difference per marker", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    m <- .cellMembership(sce, 2L, codes, c("b", "a"), rs = c(1L, 2L))
    x <- SummarizedExperiment::assay(sce)[c("b", "a"), 2]
    expect_equal(attr(m, "own"), "1")
    # Competitor is the nearest node outside the selection.
    d2 <- colSums((t(codes) - x)^2)
    expect_equal(attr(m, "competitor"), names(which.min(d2[c("3", "4")])))
    expect_equal(sum(m$contribution),
        attr(m, "d2Competitor") - attr(m, "d2Own"))
    expect_equal(m$value, unname(x))
    expect_null(.cellMembership(sce, 2L, codes, "nope", rs = 1L))
})

test_that("plotRidgeInteractive draws rows, band shapes and clickable dots", {
    sce <- make_interactive_sce()
    p <- plotRidgeInteractive(sce, c("b", "a"), rs = c(1L, 2L),
        coverage = 0.9, clusterMarkers = "a")
    b <- plotly::plotly_build(p)
    nms <- vapply(b$x$data, function(t) {
        if (is.null(t$name)) "" else t$name
    }, character(1))
    expect_true(all(c("all cells", "selection") %in% nms))
    dots <- b$x$data[grepl("^outlier", nms)]
    expect_true(length(dots) > 0L)
    expect_true(all(vapply(dots, function(t) !is.null(t$customdata),
        logical(1))))
    # One band rect per row plus two mean lines per row.
    shapes <- b$x$layout$shapes
    expect_equal(sum(vapply(shapes, function(s) s$type == "rect",
        logical(1))), 2L)
    # Sorted rows: marker "a" (larger shift) on top, marked as clustering.
    labs <- vapply(b$x$layout$annotations, function(a) a$text, character(1))
    expect_true(any(grepl("<b>a\\*</b>", labs)))
})

test_that("plotRidgeInteractive highlights the clicked cell on every row", {
    sce <- make_interactive_sce()
    p <- plotRidgeInteractive(sce, c("b", "a"), rs = c(1L, 2L),
        coverage = 0.9, highlightCell = 1L)
    b <- plotly::plotly_build(p)
    tr <- Filter(function(t) identical(t$name, "clicked cell"), b$x$data)
    expect_length(tr, 1L)
    expect_length(tr[[1]]$x, 2L)
    # On a shared axis, x is the real value of the cell per (sorted) row.
    expect_equal(sort(tr[[1]]$x),
        sort(unname(SummarizedExperiment::assay(sce)[, 1])))
})

test_that("plotRidgeInteractive scales rows separately for mixed scales", {
    sce <- make_interactive_sce()
    SummarizedExperiment::assay(sce)["b", ] <-
        SummarizedExperiment::assay(sce)["b", ] * 1e4
    b <- plotly::plotly_build(plotRidgeInteractive(sce, c("b", "a"),
        rs = c(1L, 2L), coverage = 0.9))
    xs <- unlist(lapply(b$x$data, function(t) t$x))
    expect_true(all(xs >= -1e-8 & xs <= 1 + 1e-8))
    expect_false(isTRUE(b$x$layout$xaxis$showticklabels))
})

test_that("plotRidgeInteractive returns a message plot for empty input", {
    sce <- make_interactive_sce()
    for (p in list(plotRidgeInteractive(sce, "nope", rs = 1L),
        plotRidgeInteractive(sce, "a", rs = 99L))) {
        expect_s3_class(p, "plotly")
    }
    expect_s3_class(plotCellMembership(NULL), "plotly")
})

test_that("server ridge click selects a cell and explains its membership", {
    test_app <- make_test_app()
    shiny::testServer(
        app = test_app$app,
        expr = {
            suppressWarnings(session$setInputs(clusterNumbers = "1,2"))
            session$flushReact()
            session$elapse(1100)
            suppressWarnings(session$flushReact())
            ridgeServer <- get("ridgeServer",
                envir = environment(session$testserver_serverfn)
            )
            expect_null(ridgeServer$ridgeCell())
            ridgeServer$ridgeCell(5L)
            session$flushReact()
            expect_true(!is.null(output$ridgeCellWhy))
            # A new node selection clears the clicked cell.
            rsUsed(c(3L))
            session$flushReact()
            expect_null(ridgeServer$ridgeCell())
        }
    )
})

test_that("plotCellMembership titles the own and competing node", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    m <- .cellMembership(sce, 2L, codes, c("b", "a"), rs = 1L)
    b <- plotly::plotly_build(plotCellMembership(m, "cell 2"))
    title <- b$x$layout$title$text
    expect_true(grepl("own node 1", title, fixed = TRUE))
    expect_true(grepl(attr(m, "competitor"), title, fixed = TRUE))
    bars <- b$x$data[[1]]
    expect_equal(sort(bars$x), sort(m$contribution))
    # Positive bars (pull to own node) and negative bars differ in colour.
    cols <- bars$marker$color
    expect_equal(length(unique(cols[m$contribution[match(
        as.character(bars$y), m$marker)] >= 0])), 1L)
})

test_that("plotCellMembership warns when the own node is not the nearest", {
    sce <- make_interactive_sce()
    codes <- S4Vectors::metadata(sce)$SOM_codes
    # Shift node 1's code far away: node 1 is then not the nearest code.
    codes["1", ] <- codes["1", ] + 100
    m <- .cellMembership(sce, 2L, codes, c("b", "a"), rs = 1L)
    expect_false(attr(m, "isNearest"))
    title <- plotly::plotly_build(plotCellMembership(m))$x$layout$title$text
    expect_true(grepl("is not the nearest", title, fixed = TRUE))
})

test_that("plotRidgeInteractive draws per-node lines only for picked nodes", {
    sce <- make_interactive_sce()
    nodeNames <- function(p) {
        b <- plotly::plotly_build(p)
        nms <- vapply(b$x$data, function(t) {
            if (is.null(t$name)) "" else t$name
        }, character(1))
        unique(nms[grepl("^node ", nms)])
    }
    # Two nodes (<= maxClusters): both get their own line.
    expect_setequal(nodeNames(plotRidgeInteractive(sce, c("a", "b"),
        rs = c(1L, 2L))), c("node 1", "node 2"))
    # Above the cap only the pooled selection is drawn...
    expect_length(nodeNames(plotRidgeInteractive(sce, c("a", "b"),
        rs = c(1L, 2L), maxClusters = 1L)), 0L)
    # ...unless nodes are picked explicitly.
    expect_equal(nodeNames(plotRidgeInteractive(sce, c("a", "b"),
        rs = c(1L, 2L), maxClusters = 1L, clusterPick = 2L)), "node 2")
})

test_that("plotRidgeInteractive keeps the input row order when not sorting", {
    sce <- make_interactive_sce()
    b <- plotly::plotly_build(plotRidgeInteractive(sce, c("b", "a"),
        rs = c(1L, 2L), sortRows = FALSE))
    labs <- vapply(b$x$layout$annotations, function(a) a$text, character(1))
    ys <- vapply(b$x$layout$annotations, function(a) {
        if (identical(a$yref, "y")) a$y else NA_real_
    }, numeric(1))
    rowLab <- !is.na(ys) & grepl("sd$", labs)
    # First marker ("b") on top.
    expect_true(startsWith(labs[rowLab][which.max(ys[rowLab])], "b"))
})

test_that("ridge picker offers the selected nodes, not cells indexed by id", {
    # Regression: choices were cluster_id[rs], i.e. the clusters of the
    # cells whose *column index* equals a selected node id.
    sce <- make_interactive_sce()
    sce$cluster_id <- factor(rep(c(4L, 3L, 2L, 1L, 5L, 6L, 7L, 8L),
        each = 25))
    ui <- as.character(.ridgeClusterUI(sce, rs = c(2L, 3L, 5L, 6L, 7L, 8L),
        current = NULL))
    for (id in c(2, 3, 5, 6, 7, 8)) {
        expect_true(grepl(sprintf("value=\"%d\"", id), ui, fixed = TRUE))
    }
    expect_false(grepl("value=\"4\"", ui, fixed = TRUE))
    expect_null(.ridgeClusterUI(sce, rs = c(2L, 3L), current = NULL))
})

test_that("ridge plots cap individual node lines", {
    sce <- make_interactive_sce()
    cid <- rep(seq_len(40L), each = 5L)
    sce$cluster_id <- factor(cid)
    b <- plotly::plotly_build(plotRidgeInteractive(sce, c("a", "b"),
        rs = seq_len(40L), clusterPick = seq_len(40L)))
    nms <- unique(vapply(b$x$data, function(t) {
        if (is.null(t$name)) "" else t$name
    }, character(1)))
    expect_equal(sum(grepl("^node ", nms)), .ridgeMaxLineClusters)
    expect_length(.ridgeLineClusters(cid, seq_along(cid), seq_len(40L), 5L),
        .ridgeMaxLineClusters)
})
