# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# plot-ridge-trace.R ----
# "Marker ridge" box: one row per selected marker on a shared x-axis.
# Per row: filled pastel density of all cells, solid density line for cells
# of the currently selected SOM nodes, dashed threshold ticks at
# mean +/- threshold * sd, and capped "trace" lines connecting the value of
# selected-cluster cells across the marker rows. Only cells that are extreme
# (|value - mean| > threshold * sd) on at least one marker are traced, and
# each trace is coloured by the marker on which the cell deviates most -
# this answers "why is this cell part of the selection".

#' Compute Per-Marker Densities and Statistics Over All Cells
#'
#' Builds the shared-x-axis density data used by the marker ridge box. For
#' every marker a density is estimated with \code{stats::density()} on a
#' common grid spanning all selected markers, together with the per-marker
#' mean and sd used by the threshold rule.
#'
#' @param sce A \code{SingleCellExperiment}; row \code{i} of the first assay
#'   is the expression of marker \code{i} across all cells.
#' @param somCodesName Name of the metadata slot holding the SOM codes matrix;
#'   only used to fall back to \code{rownames(sce)} when no markers are given.
#' @param markers Character vector of marker names (rownames of the assay).
#' @param nGrid Number of grid points for the density estimates.
#'
#' @return A list with:
#'   \describe{
#'     \item{allDensity}{Data frame with \code{marker}, \code{x}, \code{y}.}
#'     \item{markerStats}{Named list of \code{c(mean = , sd = )} per marker.}
#'   }
#'
#' @importFrom ggplot2 ggplot labs xlab ylab scale_y_continuous
#'     scale_fill_manual scale_colour_manual scale_linewidth_identity
#'     scale_alpha_identity geom_segment geom_text theme_void
#' @importFrom stats approx complete.cases density sd
#' @importFrom grid gpar grobName grobTree nullGrob pointsGrob polylineGrob
#' @keywords internal
.ridgeAllDensities <- function(sce,
    somCodesName = "SOM_codes",
    markers,
    nGrid = 512L) {
    mat <- SummarizedExperiment::assays(sce)[[1]]
    values <- mat[markers, , drop = FALSE]
    xRange <- range(values[!is.na(values)])
    grid <- seq(xRange[1], xRange[2], length.out = nGrid)

    allDensity <- vector("list", length(markers))
    markerStats <- vector("list", length(markers))
    for (i in seq_along(markers)) {
        v <- values[i, ]
        v <- v[!is.na(v)]
        d <- stats::density(v, n = nGrid, from = xRange[1], to = xRange[2])
        allDensity[[i]] <- data.frame(
            marker = markers[i],
            x = grid,
            y = d$y
        )
        markerStats[[i]] <- c(mean = mean(v), sd = stats::sd(v))
    }
    names(markerStats) <- markers

    allDensity <- do.call(rbind, allDensity)
    allDensity$marker <- factor(allDensity$marker, levels = markers)
    list(allDensity = allDensity, markerStats = markerStats)
}

#' Extract Trace Lines for Extreme Cells of the Selected Clusters
#'
#' For the cells belonging to the selected SOM nodes (\code{rs}) that deviate
#' more than \code{threshold} standard deviations from the mean on at least
#' one marker, builds the polyline coordinates across all selected markers
#' plus the metadata needed for colouring and capping.
#'
#' @inheritParams .ridgeAllDensities
#' @param rs Integer vector of selected SOM node ids; cells whose
#'   \code{cluster_id} is in \code{rs} are considered.
#' @param threshold Numeric; cells deviate more than \code{threshold} standard
#'   deviations count as extreme.
#' @param maxLines Integer; at most this many traces are drawn, keeping the
#'   cells with the largest maximum absolute z-score.
#'
#' @return A list with:
#'   \describe{
#'     \item{traces}{Data frame with \code{cellId}, \code{marker}, \code{x}
#'       (cell value), \code{isExtreme}; one row per cell and selected marker,
#'       sorted so that \code{geom_path} connects along the marker axis.}
#'     \item{traceMeta}{Data frame with \code{cellId}, \code{extremeMarker},
#'       \code{maxAbsZ} and \code{cluster} for the traced cells.}
#'   }
#'   Returns \code{list(traces = NULL, traceMeta = NULL)} when no cells
#'   qualify.
#'
#' @keywords internal
.ridgeTraceData <- function(sce, somCodesName = "SOM_codes", markers, rs,
    threshold = 0.99, maxLines = 50L) {
    mat <- SummarizedExperiment::assays(sce)[[1]]
    cid <- as.integer(as.character(
        SummarizedExperiment::colData(sce)$cluster_id
    ))
    sel <- which(cid %in% rs)
    ext <- .ridgeExtremeCells(mat, markers, sel, threshold, maxLines)
    if (is.null(ext)) {
        return(list(traces = NULL, traceMeta = NULL))
    }
    cellIds <- ext$cellIds
    # One row per cell and marker, sorted for geom_path connections.
    traces <- expand.grid(
        cellId = cellIds,
        marker = factor(markers, levels = markers),
        stringsAsFactors = FALSE
    )
    traces <- traces[order(traces$cellId, traces$marker), , drop = FALSE]
    # Index positionally: cell colnames can be missing or shuffled in real
    # SCEs, and name-based subsetting would silently misalign. The columns
    # of `sub` follow `sel`, so cells are matched against `sel`.
    # as.matrix() also keeps data.frame / DelayedMatrix assays working.
    sub <- as.matrix(mat[markers, sel, drop = FALSE])
    traces$x <- sub[cbind(
        match(as.character(traces$marker), markers),
        match(traces$cellId, sel)
    )]
    traces$isExtreme <- traces$marker %in% ext$extremeMarker[match(
        traces$cellId, cellIds)]
    traceMeta <- data.frame(
        cellId = cellIds,
        extremeMarker = ext$extremeMarker,
        maxAbsZ = ext$maxAbsZ,
        cluster = as.character(
            SummarizedExperiment::colData(sce)$cluster_id[cellIds]
        ),
        stringsAsFactors = FALSE
    )
    list(traces = traces, traceMeta = traceMeta)
}

#' Selected Cells Beyond the All-Cell Threshold
#'
#' @param mat Marker x cell matrix (first assay).
#' @param markers Markers to test.
#' @param sel Column indices of the selected cells.
#' @param threshold Absolute all-cell z-score above which a cell counts as
#'   extreme on a marker.
#' @param maxLines Keep at most this many cells (the most extreme).
#'
#' @return \code{NULL} when no cell qualifies, otherwise a list with
#'   \code{cellIds}, their \code{extremeMarker} (marker of maximal |z|) and
#'   \code{maxAbsZ}.
#'
#' @keywords internal
.ridgeExtremeCells <- function(mat, markers, sel, threshold, maxLines) {
    if (length(sel) == 0L) {
        return(NULL)
    }
    stats <- vapply(markers, function(m) {
        v <- mat[m, ]
        v <- v[!is.na(v)]
        c(mean = mean(v), sd = stats::sd(v))
    }, numeric(2))
    mu <- unname(stats["mean", markers, drop = TRUE])
    sds <- unname(stats["sd", markers, drop = TRUE])
    z <- abs((as.matrix(mat[markers, sel, drop = FALSE]) - mu) / sds)
    isExtreme <- apply(z, 2, function(col) any(col > threshold))
    if (!any(isExtreme)) {
        return(NULL)
    }
    zSel <- z[, isExtreme, drop = FALSE]
    cellIds <- sel[isExtreme]
    extremeMarker <- apply(zSel, 2, function(col) markers[which.max(col)])
    maxAbsZ <- apply(zSel, 2, max)
    # Cap: keep the most extreme cells.
    keep <- if (length(cellIds) > maxLines) {
        order(-maxAbsZ)[seq_len(maxLines)]
    } else {
        seq_along(cellIds)
    }
    list(cellIds = cellIds[keep], extremeMarker = unname(extremeMarker[keep]),
        maxAbsZ = unname(maxAbsZ[keep]))
}


#' Build the Marker Ridge Trace Plot
#'
#' Renders the "marker ridge" plot: one row per selected marker, the density
#' of all cells as a filled pastel area, the density of the selected-cluster
#' cells as a solid line, dashed threshold ticks at mean +/- threshold * sd,
#' and the capped trace lines of extreme cells connected across rows.
#'
#' @param sce A \code{SingleCellExperiment} with the full dataset.
#' @param somCodesName Name of the metadata slot holding the SOM codes.
#' @param markers Character vector of markers (assay rownames) to show.
#' @param rs Integer vector of selected SOM node ids.
#' @param threshold Numeric; the sd threshold (default 0.99).
#' @param maxLines Integer; cap on the number of trace lines.
#' @param maxClusters Integer; at most this many per-cluster density lines
#'   are drawn; with more selected clusters only the combined line shows.
#' @param clusterPick Integer vector of cluster ids to draw individually when
#'   more than \code{maxClusters} clusters are selected; \code{NULL} (the
#'   default) keeps the automatic cap behaviour.
#' @param showTraces Logical; when \code{FALSE} the connecting trace lines
#'   between marker rows are dropped and only the outlier dots remain
#'   (default \code{TRUE} keeps the lines).
#' @param nGrid Number of grid points for the density estimates.
#'
#' @return A \code{ggplot} object; a placeholder plot when there is nothing
#'   to show.
#'
#' @keywords internal
plotRidgeTraceFunc <- function(sce, somCodesName = "SOM_codes", markers, rs,
    threshold = 0.99, maxLines = 50L, maxClusters = 5L, clusterPick = NULL,
    nGrid = 512L, showTraces = TRUE) {
    markers <- intersect(markers, rownames(sce))
    if (length(markers) == 0L) {
        return(.ridgePlaceholder("No valid markers selected"))
    }
    ridge <- .ridgeAllDensities(sce, somCodesName, markers, nGrid = nGrid)
    if (length(rs) == 0L) {
        return(.ridgePlaceholder("No nodes selected"))
    }
    traces <- .ridgeTraceData(sce, somCodesName, markers, rs = rs,
        threshold = threshold, maxLines = maxLines)
    cid <- as.integer(as.character(
        SummarizedExperiment::colData(sce)$cluster_id
    ))
    sel <- which(cid %in% rs)
    lineClusters <- .ridgeLineClusters(cid, sel, clusterPick, maxClusters)
    clusterNames <- as.character(lineClusters)
    # Whole brewer colours, recycled: colorRampPalette(Set2)(n) interpolates
    # across the ramp and produces washed-out greys for small n.
    clusterCols <- stats::setNames(rep_len(
        RColorBrewer::brewer.pal(8L, "Set2"), length(clusterNames)
    ), clusterNames)
    markerCols <- stats::setNames(rep_len(
        RColorBrewer::brewer.pal(8L, "Set2"), length(markers)
    ), markers)
    selDensity <- .ridgeSelDensities(SummarizedExperiment::assays(sce)[[1]],
        markers, sel, cid, lineClusters, ridge, nGrid)
    geom <- .ridgeRowGeometry(ridge$allDensity, selDensity)
    areaData <- .ridgeAreaData(ridge$allDensity, selDensity, markers,
        clusterNames, geom)
    areaData <- .ridgeAppendTraces(areaData, traces, selDensity, ridge,
        markers, clusterNames, geom, showTraces)
    p <- .ridgeTraceGgplot(areaData, c(markerCols, clusterCols))
    .ridgeTraceOverlays(p, markers, ridge$markerStats, threshold,
        geom$rowHeight, c(markerCols, clusterCols))
}

#' Placeholder ggplot With a Title
#'
#' @param title Title text.
#'
#' @return An empty ggplot.
#'
#' @keywords internal
.ridgePlaceholder <- function(title) {
    ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = title)
}

#' Nodes Drawn as Individual Density Lines
#'
#' Individual lines are capped at \code{maxClusters}; with more selected
#' nodes only the combined line is drawn, unless nodes were picked
#' explicitly via \code{clusterPick}. At most \code{.ridgeMaxLineClusters}
#' lines are drawn in any case.
#'
#' @param cid Integer cluster id of every cell.
#' @param sel Column indices of the selected cells.
#' @param clusterPick Picked node ids, or \code{NULL}.
#' @param maxClusters Cap on the number of individual lines.
#'
#' @return Integer node ids (possibly empty).
#'
#' @keywords internal
.ridgeLineClusters <- function(cid, sel, clusterPick, maxClusters) {
    clusterIds <- sort(unique(cid[sel]))
    lines <- if (!is.null(clusterPick)) {
        intersect(clusterPick, clusterIds)
    } else if (length(clusterIds) <= maxClusters) {
        clusterIds
    } else {
        integer(0)
    }
    # Hard cap: hundreds of lines are unreadable and overflow the widget.
    utils::head(lines, .ridgeMaxLineClusters)
}

#' Densities of the Selection: Combined and per Node
#'
#' @param mat Marker x cell matrix (first assay).
#' @param markers Markers (rows).
#' @param sel Column indices of the selected cells.
#' @param cid Integer cluster id of every cell.
#' @param lineClusters Nodes drawn as individual lines.
#' @param ridge Output of \code{.ridgeAllDensities()} (shared grid).
#' @param nGrid Number of grid points.
#'
#' @return A data frame with \code{marker}, \code{x}, \code{y},
#'   \code{density} (= \code{y}), \code{cluster} ("combined" or the node id)
#'   and \code{linewidth}; \code{NULL} when nothing can be estimated.
#'
#' @keywords internal
.ridgeSelDensities <- function(mat, markers, sel, cid, lineClusters, ridge,
    nGrid) {
    from <- min(ridge$allDensity$x)
    to <- max(ridge$allDensity$x)
    densFrame <- function(cells, label) {
        parts <- lapply(markers, function(m) {
            vm <- as.numeric(mat[m, cells])
            vm <- vm[!is.na(vm)]
            if (length(vm) < 2L) {
                return(NULL)
            }
            data.frame(marker = m,
                x = ridge$allDensity$x[ridge$allDensity$marker == m],
                y = stats::density(vm, n = nGrid, from = from, to = to)$y)
        })
        out <- do.call(rbind, parts)
        if (!is.null(out)) out$cluster <- label
        out
    }
    selDensity <- c(
        list(if (length(sel) >= 2L) densFrame(sel, "combined")),
        lapply(lineClusters, function(k) {
            densFrame(sel[cid[sel] == k], as.character(k))
        })
    )
    selDensity <- do.call(rbind, selDensity)
    if (!is.null(selDensity)) {
        selDensity$marker <- factor(selDensity$marker, levels = markers)
        selDensity$density <- selDensity$y
        selDensity$linewidth <- ifelse(
            selDensity$cluster == "combined", 1.1, 0.7
        )
    }
    selDensity
}

#' Row Height and Cluster Stacking Step of the Ridge Rows
#'
#' Joyplot style: within each marker row, each selected node is its own
#' semi-transparent filled ridge, stacked on the row base (offset per
#' node), so overlapping distributions remain readable. The row height must
#' cover the tallest band in a row: a node density can peak higher than
#' the all-cell density, so the maximum over every curve plus the stacking
#' steps sets the scale, with 25\% headroom between rows.
#'
#' @param allDensity All-cell densities (\code{.ridgeAllDensities()}).
#' @param selDensity Output of \code{.ridgeSelDensities()}, or \code{NULL}.
#'
#' @return A list with \code{rowHeight} and \code{clusterStep}.
#'
#' @keywords internal
.ridgeRowGeometry <- function(allDensity, selDensity) {
    allPeak <- max(allDensity$y)
    clusterStep <- allPeak * 0.08
    stackTop <- 0
    if (!is.null(selDensity) && nrow(selDensity) > 0L) {
        isComb <- selDensity$cluster == "combined"
        if (any(!isComb)) {
            nStack <- length(unique(selDensity$cluster[!isComb]))
            stackTop <- max(selDensity$y[!isComb]) + nStack * clusterStep
        }
        stackTop <- max(stackTop, max(selDensity$y[isComb]))
    }
    list(rowHeight = max(allPeak, stackTop) * 1.25, clusterStep = clusterStep)
}

#' Ridge Rows of the Trace Plot in ggridges Encoding
#'
#' Background ridges (all cells) and selection ridges, offset per marker
#' row (first marker on top) and, for individual nodes, per stacking step.
#' In the ggridges identity encoding \code{y} is the ridge base and
#' \code{height} stacks upward.
#'
#' @param allDensity All-cell densities.
#' @param selDensity Output of \code{.ridgeSelDensities()}, or \code{NULL}.
#' @param markers Markers (rows).
#' @param clusterNames Node ids drawn individually, as character.
#' @param geom Output of \code{.ridgeRowGeometry()}.
#'
#' @return A data frame of ridgeline rows with \code{fillKey},
#'   \code{datatype} and one \code{ridgeGroup} per (marker, node) ridge.
#'
#' @keywords internal
.ridgeAreaData <- function(allDensity, selDensity, markers, clusterNames,
    geom) {
    rowBase <- function(m) {
        geom$rowHeight * (length(markers) - match(as.character(m), markers))
    }
    bg <- allDensity
    bg$ymin <- rowBase(bg$marker)
    bg$ymax <- bg$y + bg$ymin
    bg[c("cluster", "linewidth", "lineColour", "alpha", "linetype")] <-
        list(NA_character_, 0.5, "bgline", 0.25, 1)
    bg$y <- bg$ymin
    bg$height <- bg$ymax - bg$ymin
    keep <- c("marker", "x", "ymin", "ymax", "cluster", "linewidth",
        "linetype", "y", "height", "alpha", "lineColour")
    areaData <- bg[, keep]
    if (!is.null(selDensity) && nrow(selDensity) > 0L) {
        sd <- selDensity
        sd$cluster <- factor(sd$cluster, levels = c("combined", clusterNames))
        # The combined ridge sits at the row base; individual nodes are
        # raised one step each and stay inside their own band.
        offset <- ifelse(sd$cluster == "combined", 0,
            geom$clusterStep * match(as.character(sd$cluster), clusterNames))
        sd$ymin <- rowBase(sd$marker) + offset
        sd$ymax <- sd$ymin + sd$density
        sd$y <- sd$ymin
        sd$height <- sd$ymax - sd$ymin
        sd[c("alpha", "linetype")] <- list(0.25, 1)
        sd$lineColour <- as.character(sd$cluster)
        areaData <- rbind(areaData, sd[, keep])
    }
    areaData$fillKey <- ifelse(is.na(areaData$cluster),
        as.character(areaData$marker), as.character(areaData$cluster))
    areaData$datatype <- "ridgeline"
    ridgeKeys <- paste(as.character(areaData$marker), ifelse(
        is.na(areaData$cluster), "bg", as.character(areaData$cluster)))
    areaData$ridgeGroup <- match(ridgeKeys, unique(ridgeKeys))
    areaData
}

#' Append the Trace and Dot Rows to the Ridge Data
#'
#' Traces are kinked polylines across the marker rows, coloured by the
#' traced cell's node; each vertex sits on the ridge of the cell's node (or
#' the combined ridge) at the cell's value. They join the same
#' GeomRidgelineTraced layer as "path" rows (polylines) and "point" rows
#' (dots on the marker where the cell is most extreme); one group per
#' traced cell. \code{showTraces = FALSE} keeps only the dots.
#'
#' @param areaData Output of \code{.ridgeAreaData()}.
#' @param traces Output of \code{.ridgeTraceData()}.
#' @param selDensity Output of \code{.ridgeSelDensities()}, or \code{NULL}.
#' @param ridge Output of \code{.ridgeAllDensities()}.
#' @param markers Markers (rows).
#' @param clusterNames Node ids drawn individually, as character.
#' @param geom Output of \code{.ridgeRowGeometry()}.
#' @param showTraces Logical; draw the connecting polylines.
#'
#' @return \code{areaData} with the trace rows appended.
#'
#' @keywords internal
.ridgeAppendTraces <- function(areaData, traces, selDensity, ridge, markers,
    clusterNames, geom, showTraces) {
    tr <- traces$traces
    if (is.null(tr) || nrow(tr) == 0L) {
        return(areaData)
    }
    tr$cluster <- traces$traceMeta$cluster[
        match(tr$cellId, traces$traceMeta$cellId)]
    tr$y <- vapply(seq_len(nrow(tr)), function(r) {
        .ridgeTraceY(as.character(tr$marker[r]), tr$x[r],
            as.character(tr$cluster[r]), selDensity, ridge$allDensity,
            markers, clusterNames, geom)
    }, numeric(1))
    # ymin just below the bottom row base makes the draw_panel ordering
    # (groups sorted by ymin, descending) draw the traces on top.
    tr$ymin <- -0.05 * geom$rowHeight
    tr$ymax <- tr$y
    tr$lineColour <- as.character(tr$cluster)
    style <- list(height = 0, fillKey = "dot", linewidth = 0.9, alpha = 0.5,
        linetype = 1, point_shape = 21, point_size = 2, point_stroke = 1,
        point_fill = "black", point_alpha = 1)
    tr[names(style)] <- style
    tr$ridgeGroup <- max(areaData$ridgeGroup) +
        match(tr$cellId, unique(tr$cellId))
    for (col in setdiff(colnames(tr), colnames(areaData))) {
        areaData[[col]] <- NA
    }
    dots <- tr[tr$isExtreme, , drop = FALSE]
    dots$datatype <- "point"
    if (!isTRUE(showTraces)) {
        return(rbind(areaData, dots[, colnames(areaData)]))
    }
    tr$datatype <- "path"
    rbind(areaData, tr[, colnames(areaData)], dots[, colnames(areaData)])
}

#' Vertex Height of One Trace Point
#'
#' The density is interpolated at the cell's value: a value almost never
#' falls exactly on the density grid (exact matching flattened every trace
#' onto the row base).
#'
#' @param m Marker.
#' @param xv Cell value.
#' @param cl Node id of the traced cell.
#' @param selDensity Output of \code{.ridgeSelDensities()}, or \code{NULL}.
#' @param allDensity All-cell densities.
#' @param markers Markers (rows).
#' @param clusterNames Node ids drawn individually.
#' @param geom Output of \code{.ridgeRowGeometry()}.
#'
#' @return Plot y coordinate.
#'
#' @keywords internal
.ridgeTraceY <- function(m, xv, cl, selDensity, allDensity, markers,
    clusterNames, geom) {
    rowBase <- geom$rowHeight * (length(markers) - match(m, markers))
    if (is.na(rowBase)) {
        return(NA_real_)
    }
    ci <- match(cl, clusterNames)
    cOff <- if (!is.na(ci)) geom$clusterStep * ci else 0
    hasSel <- !is.null(selDensity) && nrow(selDensity) > 0L
    src <- if (hasSel) {
        key <- if (!is.na(ci)) cl else "combined"
        selDensity[as.character(selDensity$cluster) == key &
            as.character(selDensity$marker) == m, , drop = FALSE]
    } else {
        allDensity[as.character(allDensity$marker) == m, , drop = FALSE]
    }
    if (nrow(src) == 0L) {
        return(rowBase + cOff)
    }
    dens <- if (hasSel) src$density else src$y
    stats::approx(src$x, dens, xout = xv, rule = 2)$y + rowBase + cOff
}

#' Base ggplot of the Ridge Trace Plot
#'
#' All-cell ridges, node ridges and traces form one GeomRidgelineTraced
#' layer: ggridges' datatype mechanism keeps the row kinds (ridgeline /
#' path / point) apart inside the layer, so the z-order is set by the geom.
#'
#' @param areaData Ridge and trace rows.
#' @param cols Named fill colours of markers and nodes.
#'
#' @return A ggplot.
#'
#' @keywords internal
.ridgeTraceGgplot <- function(areaData, cols) {
    ggplot2::ggplot() +
        geom_ridgeline_traced(
            data = areaData,
            ggplot2::aes(x = x, y = y, height = height,
                fill = fillKey, colour = lineColour,
                datatype = datatype, group = ridgeGroup,
                linewidth = linewidth, alpha = alpha,
                # ggplot2 keeps only mapped columns in the build data, so
                # the styling and metadata columns ride along as aesthetics.
                marker = marker, cluster = cluster,
                cellId = cellId, isExtreme = isExtreme,
                point_shape = point_shape, point_size = point_size,
                point_stroke = point_stroke)
        ) +
        # Per-row line widths and fill alphas pass through unchanged.
        ggplot2::scale_linewidth_identity(guide = "none") +
        ggplot2::scale_alpha_identity(guide = "none") +
        # "combined" pools the selected nodes; its ridge sits at the row
        # base on top of the all-cell background, hence an explicit fill.
        ggplot2::scale_fill_manual(
            values = c(cols, combined = "grey20", dot = "black"),
            guide = "none"
        ) +
        cowplot::theme_cowplot(font_size = 12) +
        ggplot2::theme(
            legend.position = "none",
            panel.background = ggplot2::element_rect(fill = NA,
                color = "black")
        ) +
        ggplot2::xlab("Expression") +
        ggplot2::ylab(NULL) +
        ggplot2::scale_y_continuous(breaks = NULL)
}

#' Marker Labels, Colour Scale and Threshold Lines of the Trace Plot
#'
#' @param p Base ggplot from \code{.ridgeTraceGgplot()}.
#' @param markers Markers (rows).
#' @param markerStats Per-marker mean and sd (\code{.ridgeAllDensities()}).
#' @param threshold The sd threshold.
#' @param rowHeight Height of one marker row.
#' @param cols Named line colours of markers and nodes.
#'
#' @return The finished ggplot.
#'
#' @keywords internal
.ridgeTraceOverlays <- function(p, markers, markerStats, threshold,
    rowHeight, cols) {
    bases <- rowHeight * (length(markers) - seq_along(markers))
    # Labels at the base of each row, just inside the left panel edge
    # (outside it the panel clipping hides them).
    labelData <- data.frame(marker = factor(markers, levels = markers),
        x = -Inf, y = bases)
    # Dashed lines at mean +/- threshold * sd spanning each row, so the
    # outlying tails read as "outside the dashed lines".
    tickData <- do.call(rbind, lapply(seq_along(markers), function(i) {
        st <- markerStats[[markers[i]]]
        data.frame(marker = markers[i],
            x = st["mean"] + c(-1, 1) * threshold * st["sd"],
            y = bases[i], yend = bases[i] + rowHeight * 0.9)
    }))
    p +
        ggplot2::geom_text(data = labelData,
            ggplot2::aes(x = x, y = y, label = marker),
            hjust = -0.1, vjust = -0.3, inherit.aes = FALSE, size = 3.5) +
        # One colour scale: the combined line is black, the background
        # ridge line invisible.
        ggplot2::scale_colour_manual(
            values = c(cols, combined = "black", bgline = "transparent")
        ) +
        ggplot2::geom_segment(data = tickData,
            ggplot2::aes(x = x, xend = x, y = y, yend = yend,
                group = marker, colour = marker),
            linetype = "dashed", linewidth = 0.4)
}
