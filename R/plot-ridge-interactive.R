# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# plot-ridge-interactive.R ----
# Interactive (plotly) "marker ridge" box. One row per marker on a shared
# x-axis:
#   * faint grey area  - density of all cells,
#   * blue area/line   - density of the selected cells (selected SOM nodes),
#   * orange band      - central `coverage` fraction of the selected cells
#                        (empirical quantiles, e.g. 0.5% - 99.5% for 0.99),
#   * solid blue line  - selection mean; dotted grey line - all-cell mean,
#   * dots             - selected cells outside the band on that marker
#                        (the remaining 1 - coverage), hover for details,
#                        click to trace the cell across every row and to
#                        explain its SOM node membership.
# The static ggplot version (plotRidgeTraceFunc) is kept for the PDF export.

#' Per-Marker Statistics of All Cells and of the Selected Cells
#'
#' @param mat Numeric matrix, markers in rows and cells in columns.
#' @param sel Integer column indices of the selected cells.
#' @param coverage Numeric in (0, 1]; fraction of the selected cells inside
#'   the band. The band runs from the \code{(1 - coverage) / 2} to the
#'   \code{(1 + coverage) / 2} quantile of the selected cells, so about
#'   \code{1 - coverage} of them fall outside on each marker.
#'
#' @return A data frame with one row per marker (in the row order of
#'   \code{mat}): \code{marker}, \code{allMean}, \code{allSd},
#'   \code{selMean}, \code{selSd}, \code{lower}, \code{upper} (band
#'   limits) and \code{separation}, the selection mean shift in all-cell sd
#'   units.
#'
#' @keywords internal
.ridgeSelectionStats <- function(mat, sel, coverage = 0.99) {
    coverage <- min(max(coverage, 0), 1)
    probs <- c((1 - coverage) / 2, (1 + coverage) / 2)
    rowStats <- function(m) {
        vapply(seq_len(nrow(m)), function(i) {
            v <- m[i, ]
            v <- v[!is.na(v)]
            if (length(v) == 0L) {
                return(c(NA_real_, NA_real_))
            }
            c(mean(v), if (length(v) > 1L) stats::sd(v) else NA_real_)
        }, numeric(2))
    }
    all <- rowStats(mat)
    selS <- rowStats(mat[, sel, drop = FALSE])
    out <- data.frame(
        marker = rownames(mat),
        allMean = all[1, ], allSd = all[2, ],
        selMean = selS[1, ], selSd = selS[2, ],
        stringsAsFactors = FALSE
    )
    band <- vapply(seq_len(nrow(mat)), function(i) {
        v <- mat[i, sel]
        v <- v[!is.na(v)]
        if (length(v) < 2L) {
            return(c(NA_real_, NA_real_))
        }
        stats::quantile(v, probs, names = FALSE, type = 7)
    }, numeric(2))
    out$lower <- band[1, ]
    out$upper <- band[2, ]
    out$separation <- (out$selMean - out$allMean) / out$allSd
    out
}

#' Outlier Cells Within the Selection
#'
#' A selected cell is an outlier on a marker when its value lies outside
#' the band (\code{stats$lower}, \code{stats$upper}) for that marker, i.e.
#' outside the central \code{coverage} fraction of the selected cells.
#'
#' @param mat Numeric matrix, markers in rows and cells in columns.
#' @param sel Integer column indices of the selected cells.
#' @param stats Output of \code{.ridgeSelectionStats()} for \code{mat}.
#' @param maxCells Integer; at most this many outlier cells are kept, those
#'   with the largest absolute within-selection z-score.
#'
#' @return A data frame with one row per (outlier cell, marker on which it
#'   is outside the band): \code{cell} (column index), \code{marker},
#'   \code{value}, \code{z} (within-selection z-score) and \code{nOut}
#'   (number of markers the cell is outside on). Zero rows when none.
#'
#' @keywords internal
.ridgeOutliers <- function(mat, sel, stats, maxCells = 50L) {
    empty <- data.frame(
        cell = integer(0), marker = character(0), value = numeric(0),
        z = numeric(0), nOut = integer(0), stringsAsFactors = FALSE
    )
    ok <- !is.na(stats$lower) & !is.na(stats$selSd) & stats$selSd > 0
    if (length(sel) < 2L || !any(ok)) {
        return(empty)
    }
    # Index rows by name: stats may be re-ordered (rows sorted by
    # separation) relative to mat.
    vals <- mat[stats$marker[ok], sel, drop = FALSE]
    z <- (vals - stats$selMean[ok]) / stats$selSd[ok]
    isOut <- !is.na(vals) &
        (vals < stats$lower[ok] | vals > stats$upper[ok])
    nOut <- colSums(isOut)
    cand <- which(nOut > 0L)
    if (length(cand) == 0L) {
        return(empty)
    }
    maxZ <- apply(abs(z[, cand, drop = FALSE]), 2, max, na.rm = TRUE)
    cand <- cand[order(-maxZ)][seq_len(min(length(cand), maxCells))]
    idx <- which(isOut[, cand, drop = FALSE], arr.ind = TRUE)
    rowI <- idx[, 1]
    colI <- cand[idx[, 2]]
    data.frame(
        cell = sel[colI],
        marker = stats$marker[ok][rowI],
        value = vals[cbind(rowI, colI)],
        z = z[cbind(rowI, colI)],
        nOut = nOut[colI],
        stringsAsFactors = FALSE
    )
}

#' Explain the SOM Node Membership of One Cell
#'
#' FlowSOM assigns a cell to the node whose code vector is nearest
#' (squared Euclidean distance over the clustering markers). This compares
#' the cell's own node with the nearest node outside the current selection
#' and splits the difference of squared distances into per-marker
#' contributions: \code{(x - competitor)^2 - (x - own)^2}. Positive values
#' pull the cell towards its own node; they sum to
#' \code{d2Competitor - d2Own}.
#'
#' @param sce A \code{SingleCellExperiment}; the first assay is used.
#' @param cell Integer column index of the cell.
#' @param somCodes SOM codes matrix (nodes x markers); row names are the
#'   node ids used in \code{colData(sce)$cluster_id}.
#' @param clusterMarkers Markers the SOM was trained on.
#' @param rs Selected SOM node ids; the competitor is taken outside them
#'   (or, when every node is selected, any node other than the own one).
#'
#' @return A data frame with \code{marker}, \code{value}, \code{ownCode},
#'   \code{compCode} and \code{contribution}, carrying the attributes
#'   \code{own}, \code{competitor}, \code{d2Own}, \code{d2Competitor},
#'   \code{nearest} and \code{isNearest}; \code{NULL} when no clustering
#'   marker is available.
#'
#' @keywords internal
.cellMembership <- function(sce, cell, somCodes, clusterMarkers, rs) {
    if (is.null(rownames(somCodes))) {
        rownames(somCodes) <- as.character(seq_len(nrow(somCodes)))
    }
    clusterMarkers <- intersect(
        intersect(clusterMarkers, colnames(somCodes)), rownames(sce)
    )
    if (length(clusterMarkers) == 0L) {
        return(NULL)
    }
    codes <- somCodes[, clusterMarkers, drop = FALSE]
    x <- as.numeric(as.matrix(
        SummarizedExperiment::assays(sce)[[1]][clusterMarkers, cell,
            drop = FALSE]
    ))
    own <- as.character(SummarizedExperiment::colData(sce)$cluster_id[cell])
    d2 <- colSums((t(codes) - x)^2)
    names(d2) <- rownames(codes)
    others <- setdiff(rownames(codes), c(as.character(rs), own))
    if (length(others) == 0L) {
        others <- setdiff(rownames(codes), own)
    }
    if (length(others) == 0L || !own %in% rownames(codes)) {
        return(NULL)
    }
    comp <- others[which.min(d2[others])]
    ownCode <- as.numeric(codes[own, ])
    compCode <- as.numeric(codes[comp, ])
    out <- data.frame(
        marker = clusterMarkers,
        value = x,
        ownCode = ownCode,
        compCode = compCode,
        contribution = (x - compCode)^2 - (x - ownCode)^2,
        stringsAsFactors = FALSE
    )
    nearest <- names(d2)[which.min(d2)]
    attr(out, "own") <- own
    attr(out, "competitor") <- comp
    attr(out, "d2Own") <- unname(d2[own])
    attr(out, "d2Competitor") <- unname(d2[comp])
    attr(out, "nearest") <- nearest
    attr(out, "isNearest") <- identical(nearest, own)
    out
}

#' Interactive Marker Ridge Plot
#'
#' Plotly version of the marker ridge box; see the file header of
#' \code{plot-ridge-interactive.R} for the visual encoding. Outlier dots
#' carry the cell's column index as \code{customdata}, so a
#' \code{plotly_click} event on \code{source} identifies the cell.
#'
#' @param sce A \code{SingleCellExperiment}; the first assay is used.
#' @param markers Character vector of markers (assay row names) to show.
#' @param rs Integer vector of selected SOM node ids.
#' @param coverage Numeric in (0, 1]; fraction of the selected cells inside
#'   the band; see \code{.ridgeSelectionStats()}.
#' @param maxCells Integer; cap on the number of outlier cells shown.
#' @param maxClusters Integer; with more selected nodes only the pooled
#'   selection density is drawn, unless \code{clusterPick} is given.
#' @param clusterPick Integer node ids to draw individually, or
#'   \code{NULL} for the automatic cap.
#' @param highlightCell Column index of a cell to trace across all rows,
#'   or \code{NULL}.
#' @param clusterMarkers Markers used for SOM clustering; their row labels
#'   are marked with an asterisk.
#' @param sortRows Logical; order rows by decreasing absolute separation of
#'   the selection from all cells.
#' @param source Plotly event source id.
#' @param maxSpanRatio Numeric; when the widest marker range exceeds the
#'   narrowest by more than this factor (e.g. transformed and untransformed
#'   channels mixed), each row is drawn on its own min-to-max scale instead
#'   of one shared expression axis.
#' @param nGrid Number of grid points for the density estimates.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @keywords internal
plotRidgeInteractive <- function(sce, markers, rs, coverage = 0.99,
    maxCells = 50L, maxClusters = 5L, clusterPick = NULL,
    highlightCell = NULL, clusterMarkers = NULL, sortRows = TRUE,
    source = "ridgePlot", maxSpanRatio = 10, nGrid = 512L) {
    markers <- intersect(markers, rownames(sce))
    if (length(markers) == 0L) {
        return(.ridgeEmptyPlot("No valid markers selected", source))
    }
    cid <- as.integer(as.character(
        SummarizedExperiment::colData(sce)$cluster_id
    ))
    sel <- which(cid %in% rs)
    if (length(sel) == 0L) {
        return(.ridgeEmptyPlot("No cells in the selected nodes", source))
    }
    ctx <- .ridgeContext(sce, markers, cid, sel, coverage, sortRows,
        maxSpanRatio, nGrid, clusterPick, maxClusters, clusterMarkers)
    p <- plotly::plot_ly(source = source)
    shapes <- list()
    annotations <- list()
    for (r in seq_along(ctx$markers)) {
        p <- .ridgeRowTraces(p, ctx, r)
        deco <- .ridgeRowDecorations(ctx, r)
        shapes <- c(shapes, deco$shapes)
        annotations <- c(annotations, deco$annotations)
    }
    p <- .ridgeAddOutliers(p, ctx,
        .ridgeOutliers(ctx$mat, sel, ctx$st, maxCells))
    p <- .ridgeAddHighlight(p, ctx, highlightCell)
    .ridgeLayout(p, ctx, shapes, annotations)
}

#' Placeholder Plotly Widget With a Message
#'
#' @param msg Message to show.
#' @param source Plotly event source id.
#'
#' @return A plotly widget.
#'
#' @keywords internal
.ridgeEmptyPlot <- function(msg, source = "ridgePlot") {
    plotly::plot_ly(source = source, type = "scatter", mode = "markers") %>%
        plotly::layout(
            xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
            annotations = list(list(text = msg, showarrow = FALSE,
                xref = "paper", yref = "paper", x = 0.5, y = 0.5))
        )
}

#' Shared Data of the Interactive Ridge Plot
#'
#' Collects everything the per-row builders need: the marker matrix and
#' statistics (rows optionally sorted by separation), the row positions, the
#' per-row x ranges and the value-to-plot mapping \code{tx()}, and the
#' densities of all cells, the selection and the per-node lines per row.
#'
#' @inheritParams plotRidgeInteractive
#' @param cid Integer cluster id of every cell.
#' @param sel Column indices of the selected cells.
#'
#' @return A list (the ridge context).
#'
#' @keywords internal
.ridgeContext <- function(sce, markers, cid, sel, coverage, sortRows,
    maxSpanRatio, nGrid, clusterPick, maxClusters, clusterMarkers) {
    mat <- as.matrix(SummarizedExperiment::assays(sce)[[1]][markers, ,
        drop = FALSE])
    st <- .ridgeSelectionStats(mat, sel, coverage)
    if (isTRUE(sortRows)) {
        st <- st[order(-abs(st$separation), na.last = TRUE), , drop = FALSE]
        markers <- st$marker
    }
    ctx <- c(
        list(mat = mat, st = st, markers = markers, sel = sel, cid = cid,
            coverage = coverage, nGrid = nGrid,
            clusterMarkers = clusterMarkers,
            base = stats::setNames(length(markers) - seq_along(markers),
                markers),
            cellLabel = .ridgeCellLabeller(sce, cid)),
        .ridgeXScale(mat, markers, sel, maxSpanRatio)
    )
    clusterIds <- sort(unique(cid[sel]))
    ctx$lineClusters <- if (!is.null(clusterPick)) {
        intersect(clusterPick, clusterIds)
    } else if (length(clusterIds) <= maxClusters &&
        length(clusterIds) > 1L) {
        clusterIds
    } else {
        integer(0)
    }
    ctx$lineClusters <- utils::head(ctx$lineClusters, .ridgeMaxLineClusters)
    ctx$clusterCols <- rep_len(RColorBrewer::brewer.pal(8L, "Set2"),
        length(ctx$lineClusters))
    ctx$rows <- lapply(markers, function(m) .ridgeRowDensities(ctx, m))
    names(ctx$rows) <- markers
    ctx
}

#' X Ranges of the Ridge Rows
#'
#' Per marker the robust (0.1\% - 99.9\%) range of all cells, widened to
#' keep every selected cell visible. When the widest range exceeds the
#' narrowest by more than \code{maxSpanRatio}, each row keeps its own range
#' and is mapped onto 0 - 1; otherwise all rows share one range.
#'
#' @param mat Marker x cell matrix.
#' @param markers Marker (row) names.
#' @param sel Column indices of the selected cells.
#' @param maxSpanRatio See \code{plotRidgeInteractive()}.
#'
#' @return A list with \code{rowRange} (named list of ranges),
#'   \code{perMarker} (logical) and \code{tx}, a function
#'   \code{(marker, values)} mapping values to plot x coordinates.
#'
#' @keywords internal
.ridgeXScale <- function(mat, markers, sel, maxSpanRatio) {
    rowRange <- lapply(markers, function(m) {
        v <- mat[m, ]
        v <- v[!is.na(v)]
        r <- range(c(stats::quantile(v, c(0.001, 0.999), names = FALSE),
            range(mat[m, sel], na.rm = TRUE)))
        if (diff(r) <= 0) r + c(-0.5, 0.5) else r
    })
    names(rowRange) <- markers
    spans <- vapply(rowRange, diff, numeric(1))
    perMarker <- max(spans) / min(spans) > maxSpanRatio
    if (!perMarker) {
        rowRange[] <- list(range(unlist(rowRange)))
    }
    tx <- function(m, v) {
        if (perMarker) (v - rowRange[[m]][1]) / diff(rowRange[[m]]) else v
    }
    list(rowRange = rowRange, perMarker = perMarker, tx = tx)
}

#' Densities of One Ridge Row
#'
#' @param ctx Ridge context (see \code{.ridgeContext()}).
#' @param m Marker name.
#'
#' @return A list with the density \code{grid}, its plot coordinates
#'   \code{gx}, the densities \code{dAll}, \code{dSel}, \code{dCl} (per
#'   node line) and the common row scale \code{s}.
#'
#' @keywords internal
.ridgeRowDensities <- function(ctx, m) {
    rr <- ctx$rowRange[[m]]
    dens <- function(v) {
        v <- v[!is.na(v)]
        if (length(v) < 2L) {
            return(NULL)
        }
        stats::density(v, n = ctx$nGrid, from = rr[1], to = rr[2])$y
    }
    grid <- seq(rr[1], rr[2], length.out = ctx$nGrid)
    sel <- ctx$sel
    row <- list(
        grid = grid, gx = ctx$tx(m, grid),
        dAll = dens(ctx$mat[m, ]), dSel = dens(ctx$mat[m, sel]),
        dCl = lapply(ctx$lineClusters, function(k) {
            dens(ctx$mat[m, sel[ctx$cid[sel] == k]])
        })
    )
    # Heights are scaled per row: all curves of a row share one scale, so
    # all-cell and selection densities remain comparable within the row.
    row$s <- max(c(row$dAll, row$dSel, unlist(row$dCl), 1e-12))
    row
}

#' Hover Label Factory for Single Cells
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param cid Integer cluster id of every cell.
#'
#' @return A function mapping a column index to a label with the cell id,
#'   its node and its sample.
#'
#' @keywords internal
.ridgeCellLabeller <- function(sce, cid) {
    cd <- SummarizedExperiment::colData(sce)
    function(cell) {
        lab <- if (is.null(cd$cellId)) cell else cd$cellId[cell]
        paste0("cell ", lab, " \u00b7 node ", cid[cell],
            if (!is.null(cd$sample_id)) paste0(" \u00b7 ", cd$sample_id[cell]))
    }
}

#' Number Formatting of the Ridge Hover Texts
#'
#' @param v Numeric vector.
#'
#' @return Character vector with three significant digits.
#'
#' @keywords internal
.ridgeFmt <- function(v) formatC(v, digits = 3, format = "g")

#' Add the Density Traces of One Ridge Row
#'
#' @param p Plotly widget.
#' @param ctx Ridge context.
#' @param r Row index.
#'
#' @return The widget with the all-cell, selection and per-node traces.
#'
#' @keywords internal
.ridgeRowTraces <- function(p, ctx, r) {
    m <- ctx$markers[r]
    b <- ctx$base[[m]]
    row <- ctx$rows[[m]]
    gx <- row$gx
    ridge <- function(d) c(b + 0.85 * d / row$s, rep(b, ctx$nGrid))
    first <- r == 1L
    if (!is.null(row$dAll)) {
        p <- p %>% plotly::add_trace(type = "scatter", mode = "lines",
            x = c(gx, rev(gx)), y = ridge(row$dAll),
            fill = "toself", fillcolor = "rgba(150,150,150,0.25)",
            line = list(width = 0.5, color = "rgba(120,120,120,0.6)"),
            name = "all cells", legendgroup = "all",
            showlegend = first, hoverinfo = "skip")
    }
    if (!is.null(row$dSel)) {
        p <- p %>% plotly::add_trace(type = "scatter", mode = "lines",
            x = c(gx, rev(gx)), y = ridge(row$dSel),
            fill = "toself", fillcolor = "rgba(31,119,180,0.25)",
            line = list(width = 1.5, color = "rgb(31,119,180)"),
            name = "selection", legendgroup = "sel", showlegend = first,
            text = .ridgeRowSummary(ctx$st[r, ]), hoverinfo = "text",
            hoveron = "fills")
    }
    for (k in seq_along(ctx$lineClusters)) {
        if (is.null(row$dCl[[k]])) next
        p <- p %>% plotly::add_trace(type = "scatter", mode = "lines",
            x = gx, y = b + 0.85 * row$dCl[[k]] / row$s,
            line = list(width = 1, color = ctx$clusterCols[k]),
            name = paste("node", ctx$lineClusters[k]),
            legendgroup = paste0("node", ctx$lineClusters[k]),
            showlegend = first, hoverinfo = "skip")
    }
    p
}

#' Hover Summary of One Ridge Row
#'
#' @param st One row of \code{.ridgeSelectionStats()}.
#'
#' @return A character string.
#'
#' @keywords internal
.ridgeRowSummary <- function(st) {
    fmt <- .ridgeFmt
    paste0(
        "<b>", st$marker, "</b><br>selection mean ", fmt(st$selMean),
        " (sd ", fmt(st$selSd), ")<br>band [", fmt(st$lower), ", ",
        fmt(st$upper), "]<br>all cells mean ", fmt(st$allMean),
        " (sd ", fmt(st$allSd), ")<br>shift ",
        sprintf("%+.2f", st$separation), " sd of all cells"
    )
}

#' Band, Mean Lines and Labels of One Ridge Row
#'
#' @param ctx Ridge context.
#' @param r Row index.
#'
#' @return A list of plotly \code{shapes} and \code{annotations}.
#'
#' @keywords internal
.ridgeRowDecorations <- function(ctx, r) {
    m <- ctx$markers[r]
    b <- ctx$base[[m]]
    st <- ctx$st[r, ]
    vline <- function(x, line) {
        list(type = "line", x0 = ctx$tx(m, x), x1 = ctx$tx(m, x),
            y0 = b, y1 = b + 0.9, line = line)
    }
    shapes <- list(
        if (!is.na(st$selSd)) {
            list(type = "rect", layer = "below", line = list(width = 0),
                x0 = ctx$tx(m, st$lower), x1 = ctx$tx(m, st$upper),
                y0 = b, y1 = b + 0.9, fillcolor = "rgba(255,165,0,0.18)")
        },
        vline(st$selMean, list(color = "rgb(31,119,180)", width = 1.5)),
        vline(st$allMean, list(color = "grey", width = 1, dash = "dot"))
    )
    isClust <- m %in% ctx$clusterMarkers
    annotations <- list(list(
        text = paste0(if (isClust) "<b>", m, if (isClust) "*</b>", "<br>",
            sprintf("%+.1f sd", st$separation)),
        x = 0, xref = "paper", xanchor = "right", y = b + 0.35,
        yref = "y", showarrow = FALSE, align = "right",
        font = list(size = 11)
    ))
    if (ctx$perMarker) {
        # Real min / max of the row's own scale at its two ends.
        annotations <- c(annotations, lapply(seq_len(2L), function(e) {
            list(text = .ridgeFmt(ctx$rowRange[[m]][e]), x = e - 1,
                xref = "x", xanchor = if (e == 1L) "left" else "right",
                y = b, yref = "y", yanchor = "top", showarrow = FALSE,
                font = list(size = 9, color = "grey"))
        }))
    }
    list(shapes = Filter(Negate(is.null), shapes), annotations = annotations)
}

#' Height of the Selection Density at Given Values
#'
#' @param ctx Ridge context.
#' @param m Marker name.
#' @param v Values of marker \code{m}.
#'
#' @return Plot y coordinates on the selection density curve of the row
#'   (the row base when the row has no selection density).
#'
#' @keywords internal
.ridgeOnCurve <- function(ctx, m, v) {
    row <- ctx$rows[[m]]
    b <- ctx$base[[m]]
    if (is.null(row$dSel)) {
        return(rep(b, length(v)))
    }
    stats::approx(row$grid, b + 0.85 * row$dSel / row$s, xout = v,
        rule = 2)$y
}

#' Add the Outlier Dots to the Ridge Plot
#'
#' @param p Plotly widget.
#' @param ctx Ridge context.
#' @param outl Output of \code{.ridgeOutliers()}.
#'
#' @return The widget with one trace of dots above and one below the band.
#'
#' @keywords internal
.ridgeAddOutliers <- function(p, ctx, outl) {
    if (nrow(outl) == 0L) {
        return(p)
    }
    pos <- function(f) {
        vapply(seq_len(nrow(outl)), function(i) {
            f(outl$marker[i], outl$value[i])
        }, numeric(1))
    }
    outl$x <- pos(ctx$tx)
    outl$y <- pos(function(m, v) .ridgeOnCurve(ctx, m, v))
    outl$txt <- paste0(
        vapply(outl$cell, ctx$cellLabel, character(1)), "<br><b>",
        outl$marker, "</b> = ", .ridgeFmt(outl$value),
        "<br>z = ", sprintf("%+.2f", outl$z),
        " (selection sd)<br>outside the band on ", outl$nOut,
        " marker(s)<br><i>click to trace and explain</i>"
    )
    for (high in c(TRUE, FALSE)) {
        d <- outl[if (high) outl$z > 0 else outl$z < 0, , drop = FALSE]
        if (nrow(d) == 0L) next
        p <- p %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            x = d$x, y = d$y, customdata = d$cell,
            text = d$txt, hoverinfo = "text",
            name = paste("outlier,", if (high) "above band" else "below band"),
            marker = list(size = 8,
                color = if (high) "rgb(214,39,40)" else "rgb(44,160,44)",
                line = list(color = "white", width = 1))
        )
    }
    p
}

#' Trace the Clicked Cell Across All Ridge Rows
#'
#' @param p Plotly widget.
#' @param ctx Ridge context.
#' @param cell Column index of the cell, or \code{NULL}.
#'
#' @return The widget, with a dashed black trace when \code{cell} is valid.
#'
#' @keywords internal
.ridgeAddHighlight <- function(p, ctx, cell) {
    if (is.null(cell) || !cell %in% seq_len(ncol(ctx$mat))) {
        return(p)
    }
    markers <- ctx$markers
    hv <- ctx$mat[markers, cell]
    at <- function(f) {
        vapply(seq_along(markers), function(i) f(markers[i], hv[i]),
            numeric(1))
    }
    hz <- (hv - ctx$st$selMean) / ctx$st$selSd
    p %>% plotly::add_trace(
        type = "scatter", mode = "lines+markers",
        x = at(ctx$tx), y = at(function(m, v) .ridgeOnCurve(ctx, m, v)),
        customdata = rep(cell, length(markers)), name = "clicked cell",
        text = paste0(ctx$cellLabel(cell), "<br><b>", markers,
            "</b> = ", .ridgeFmt(hv), "<br>z = ", sprintf("%+.2f", hz),
            " (selection sd)"),
        hoverinfo = "text",
        line = list(color = "black", width = 1.5, dash = "dash"),
        marker = list(symbol = "diamond", size = 10, color = "black")
    )
}

#' Final Layout of the Interactive Ridge Plot
#'
#' @param p Plotly widget.
#' @param ctx Ridge context.
#' @param shapes,annotations Plotly shapes and annotations of the rows.
#'
#' @return The finished widget.
#'
#' @keywords internal
.ridgeLayout <- function(p, ctx, shapes, annotations) {
    coverage <- ctx$coverage
    legendNote <- list(
        text = paste0(
            "band: central ", 100 * coverage, "% of the selected cells",
            " (dots: the ", 100 * (1 - coverage), "% outside)",
            " \u00b7 solid: selection mean \u00b7 dotted: all-cell mean",
            if (length(ctx$clusterMarkers)) " \u00b7 * clustering marker"
        ),
        x = 0, xref = "paper", xanchor = "left", y = 1.0, yref = "paper",
        yanchor = "bottom", showarrow = FALSE, font = list(size = 11)
    )
    xaxis <- if (ctx$perMarker) {
        list(title = paste("Expression, each row scaled to its own",
            "range (min \u2192 max, grey labels)"),
            showticklabels = FALSE, zeroline = FALSE)
    } else {
        list(title = "Expression", zeroline = FALSE)
    }
    nRows <- length(ctx$markers)
    p %>% plotly::layout(
        shapes = shapes, annotations = c(annotations, list(legendNote)),
        xaxis = xaxis,
        yaxis = list(visible = FALSE, range = c(-0.05, nRows + 0.05),
            fixedrange = TRUE),
        margin = list(l = 110, t = 30),
        legend = list(orientation = "h", y = -0.12),
        hovermode = "closest"
    )
}

#' Plot the Node Membership Explanation of One Cell
#'
#' Horizontal bars of the per-marker contributions returned by
#' \code{.cellMembership()}: positive bars pull the cell towards its own
#' SOM node, negative bars towards the competing node.
#'
#' @param memb Output of \code{.cellMembership()}, or \code{NULL}.
#' @param cellLabel Character label of the cell used in the title.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @keywords internal
plotCellMembership <- function(memb, cellLabel = "") {
    if (is.null(memb) || nrow(memb) == 0L) {
        return(plotly::plot_ly(type = "scatter", mode = "markers") %>%
            plotly::layout(
                xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                annotations = list(list(
                    text = paste("Click an outlier dot in the ridge plot to",
                        "explain why the cell belongs to its SOM node"),
                    showarrow = FALSE, xref = "paper", yref = "paper",
                    x = 0.5, y = 0.5))
            ))
    }
    own <- attr(memb, "own")
    comp <- attr(memb, "competitor")
    fmt <- function(v) formatC(v, digits = 3, format = "g")
    memb <- memb[order(memb$contribution), , drop = FALSE]
    title <- paste0(
        cellLabel, ": own node ", own, " (d\u00b2 ", fmt(attr(memb, "d2Own")),
        ") vs nearest node outside the selection ", comp, " (d\u00b2 ",
        fmt(attr(memb, "d2Competitor")), ")",
        if (isTRUE(attr(memb, "isNearest"))) "" else paste0(
            "<br><span style='color:red'>node ", own, " is not the nearest ",
            "code (", attr(memb, "nearest"), "): SOM codes and assay ",
            "may be on different scales</span>")
    )
    plotly::plot_ly(
        x = memb$contribution,
        y = factor(memb$marker, levels = memb$marker),
        type = "bar", orientation = "h",
        marker = list(color = ifelse(memb$contribution >= 0,
            "rgb(31,119,180)", "rgb(214,39,40)")),
        text = paste0("<b>", memb$marker, "</b><br>cell ", fmt(memb$value),
            "<br>node ", own, " code ", fmt(memb$ownCode),
            "<br>node ", comp, " code ", fmt(memb$compCode),
            "<br>contribution ", sprintf("%+.3g", memb$contribution)),
        hoverinfo = "text", textposition = "none"
    ) %>% plotly::layout(
        title = list(text = title, font = list(size = 12), x = 0,
            xanchor = "left"),
        xaxis = list(title = paste0(
            "\u2190 pulls towards node ", comp, "   |   pulls towards node ",
            own, " \u2192"), zeroline = TRUE),
        yaxis = list(title = ""),
        margin = list(t = 60)
    )
}
