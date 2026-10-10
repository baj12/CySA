# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# plot-helpers.R ----
# Plot building helpers extracted from the clusterSelector() server closure.
# These functions produce static ggplot2 objects that are used both in the app
# and can be tested independently.

#' Build a Base SOM Raster Heatmap
#'
#' Creates a faceted ggplot2 raster heatmap from \code{somRasterData}. The
#' resulting object is reused as the background for selected-node overlays.
#'
#' @param somRasterData Data frame with \code{x}, \code{y}, and marker columns.
#' @param colsUsed Marker columns to include.
#'
#' @return A \code{ggplot} object or \code{NULL} if \code{somRasterData} is
#'   \code{NULL}.
#'
#' @keywords internal
.buildBaseRasterGgplot <- function(somRasterData, colsUsed) {
    if (is.null(somRasterData)) {
        return(NULL)
    }

    somRasterData <- somRasterData[, c("x", "y", colsUsed), drop = FALSE]
    raster_long <- tidyr::pivot_longer(
        somRasterData,
        cols = -c("x", "y"),
        names_to = "marker",
        values_to = "value"
    )
    marker_cols <- setdiff(names(somRasterData), c("x", "y"))
    raster_long$marker <- factor(raster_long$marker, levels = marker_cols)

    ggplot(raster_long, aes(x = x, y = y, fill = value)) +
        geom_raster() +
        facet_wrap(~marker) +
        scale_fill_viridis() +
        coord_fixed() +
        theme_minimal()
}


#' Build a FlowSOM Marker Pie Plot
#'
#' @param somCodes SOM codes matrix.
#' @param rs Selected SOM node ids.
#' @param colsUsed Marker names to include.
#' @param ncol Number of columns passed to \code{facet_wrap}. \code{NULL}
#'   lets ggplot2 choose automatically.
#' @param maxPies Maximum number of pies to plot. When more nodes are selected,
#'   only the first \code{maxPies} are shown. Default is 5.
#'
#' @importFrom tidyselect all_of
#' @return A \code{ggplot} object or \code{NULL}.
#' @keywords internal
.buildFlowSOMPiePlot <- function(somCodes, rs, colsUsed, ncol = NULL,
    maxPies = 5L) {
    rs <- as.integer(rs)
    rs <- rs[!is.na(rs)]
    if (length(rs) < 1L || any(rs <= 0L) || any(rs > nrow(somCodes))) {
        return(NULL)
    }

    ## Limit to maxPies
    if (length(rs) > maxPies) {
        rs <- rs[seq_len(maxPies)]
    }

    markers <- intersect(colsUsed, colnames(somCodes))
    if (length(markers) < 1L) {
        return(NULL)
    }

    ## Validate ncol
    if (!is.null(ncol)) {
        ncol <- .asIntegerQuiet(ncol)[1L]
        if (is.na(ncol) || ncol < 1L) ncol <- NULL
    }

    long <- .pieData(somCodes, rs, markers)
    marker_cols <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8L, "Set2")
    )(length(markers))

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = "", y = expression, fill = marker)
    ) +
        ggplot2::geom_bar(stat = "identity", width = 1L) +
        ggplot2::coord_polar("y") +
        ggplot2::facet_wrap(~cluster_id, ncol = ncol) + # NULL ncol -> auto
        ggplot2::scale_fill_manual(values = marker_cols) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text    = ggplot2::element_blank(),
            axis.title   = ggplot2::element_blank(),
            panel.grid   = ggplot2::element_blank()
        )
}


#' Long-Format Marker Codes of the Pie Nodes
#'
#' @param somCodes SOM codes matrix.
#' @param rs Node ids (one pie each, in this order).
#' @param markers Markers (pie slices).
#'
#' @return A data frame with \code{cluster_id}, \code{marker} and
#'   \code{expression}.
#'
#' @keywords internal
.pieData <- function(somCodes, rs, markers) {
    expr <- as.data.frame(somCodes[rs, markers, drop = FALSE])
    expr$cluster_id <- factor(as.character(rs),
        levels = unique(as.character(rs)))
    long <- tidyr::pivot_longer(expr, cols = tidyselect::all_of(markers),
        names_to = "marker", values_to = "expression")
    long$marker <- factor(long$marker, levels = markers)
    long
}


#' Build a SOM Raster Selection Plot
#'
#' Creates a plotly-friendly ggplot showing the full SOM grid with selected
#' nodes highlighted.
#'
#' @param somRasterData Data frame containing \code{x} and \code{y} grid
#'   columns.
#' @param rs Selected SOM node ids.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
.buildSOMRasterSelectPlot <- function(somRasterData, rs) {
    # Guard: ensure required columns exist before subsetting
    required <- c("x", "y", "id")
    missing_cols <- setdiff(required, colnames(somRasterData))
    if (length(missing_cols) > 0L) {
        warning(
            ".buildSOMRasterSelectPlot: somRasterData is missing required ",
            "column(s): ",
            paste(sQuote(missing_cols), collapse = ", "),
            ". Returning NULL."
        )
        return(NULL)
    }

    data.points <- somRasterData[, required, drop = FALSE]

    rs_valid <- as.integer(unlist(rs))
    rs_valid <- rs_valid[!is.na(rs_valid) & rs_valid %in% data.points$id]

    selected <- if (length(rs_valid) > 0L) {
        data.points[data.points$id %in% rs_valid, , drop = FALSE]
    } else {
        data.points[0L, , drop = FALSE]
    }

    ggplot2::ggplot(
        data.points,
        ggplot2::aes(x = x, y = y, key = id)
    ) +
        ggplot2::geom_point() +
        ggplot2::geom_point(
            data = selected,
            ggplot2::aes(x = x, y = y),
            color = "red",
            size = 0.9,
            inherit.aes = FALSE
        )
}

#' Compute Percentile Axis Limits for Scatter Plots
#'
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param pctl Quantile between 0.5 and 1.
#'
#' @return A named list with \code{xlim} and \code{ylim}.
#'
#' @keywords internal
.computeScatterLimits <- function(x, y, pctl) {
    if (is.null(pctl)) pctl <- 0.99
    tailP <- (1 - pctl) / 2
    list(
        xlim = stats::quantile(x, probs = c(tailP, 1 - tailP), na.rm = TRUE),
        ylim = stats::quantile(y, probs = c(tailP, 1 - tailP), na.rm = TRUE)
    )
}
