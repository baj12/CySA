# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server-ridge.R ----
# Server logic of the interactive "marker ridge" box: the plotly ridge plot,
# the click observer that picks an outlier cell, and the membership
# explanation of the clicked cell.

#' Build the Marker Ridge Box Server Logic
#'
#' Registers \code{output$ridgePlot} (interactive ridge plot),
#' \code{output$ridgeCellWhy} (SOM node membership of the clicked cell)
#' and the \code{plotly_click} observer on the ridge plot. The clicked cell
#' is reset whenever the node selection changes.
#'
#' @param input,output Shiny \code{input} / \code{output}.
#' @param sce Full \code{SingleCellExperiment}.
#' @param somCodes SOM codes matrix (nodes x markers).
#' @param colsUsed Markers used for SOM clustering (default ridge markers).
#' @param rsUsed Reactive value holding the selected SOM node ids.
#' @param selectedUpdate,refresh Reactives whose invalidation re-renders
#'   the plot (group edits).
#' @param verbose Logical; passed to \code{.safeEventData()}.
#'
#' @return A list with the reactives \code{ridgePlot} (debounced plotly
#'   widget) and \code{ridgeCell} (clicked cell index or \code{NULL}),
#'   invisibly.
#'
#' @keywords internal
.buildRidgeServer <- function(input, output, sce, somCodes, colsUsed,
    rsUsed, selectedUpdate = function() NULL, refresh = function() NULL,
    verbose = FALSE) {
    ridgeCell <- shiny::reactiveVal(NULL)

    shiny::observeEvent(rsUsed(), ridgeCell(NULL), ignoreInit = TRUE)

    shiny::observeEvent(
        .safeEventData("plotly_click", source = "ridgePlot",
            verbose = verbose),
        {
            d <- .safeEventData("plotly_click", source = "ridgePlot",
                verbose = verbose)
            click <- .parseRidgeClick(d)
            if (is.null(click)) {
                return(NULL)
            }
            ridgeCell(click$cell)
        }
    )

    ridgePlot <- shiny::reactive({
        selectedUpdate()
        refresh()
        rs <- rsUsed()
        shiny::req(rs)
        .ridgePlotFromInput(sce, input, as.integer(unlist(rs)), colsUsed,
            ridgeCell())
    }) %>% shiny::debounce(500)

    output$ridgePlot <- plotly::renderPlotly({
        ridgePlot() %>% plotly::event_register("plotly_click")
    })

    output$ridgeCellWhy <- plotly::renderPlotly(
        .ridgeCellWhyPlot(sce, somCodes, colsUsed, ridgeCell(), rsUsed())
    )

    invisible(list(ridgePlot = ridgePlot, ridgeCell = ridgeCell))
}

#' Parse a Ridge Plot Click Payload
#'
#' Outlier dots carry \code{"<marker>\x1f<cell>"} in \code{customdata} so a
#' click names both the cell and the ridge row it came from. The "clicked
#' cell" overlay trace carries a bare cell index (no marker). Payloads
#' shaped as a data frame (the usual plotly event), a bare atomic value or
#' a point-object list are all accepted; anything else parses to
#' \code{NULL}.
#'
#' @param d Event data from \code{.safeEventData("plotly_click", ...)},
#'   which may be a data frame, a point-object list or an atomic value, or
#'   \code{NULL}.
#'
#' @return A list with \code{marker} (\code{NULL} when absent) and
#'   \code{cell} (integer column index into \code{sce}), or \code{NULL}.
#'
#' @keywords internal
.parseRidgeClick <- function(d) {
    if (is.null(d) || length(d) == 0L) {
        return(NULL)
    }
    if (is.data.frame(d)) {
        cd <- d[["customdata"]][1L]
    } else if (is.list(d)) {
        # A single point object or an array of point objects.
        point <- if (!is.null(d[["customdata"]])) d else d[[1L]]
        if (!is.list(point)) {
            return(NULL)
        }
        cd <- point[["customdata"]]
    } else {
        cd <- d[1L]
    }
    if (is.null(cd)) {
        return(NULL)
    }
    if (is.list(cd)) {
        cd <- cd[[1L]]
    }
    cd <- as.character(cd)
    if (length(cd) != 1L || is.na(cd) || !nzchar(cd)) {
        return(NULL)
    }
    parts <- strsplit(cd, "\x1f", fixed = TRUE)[[1]]
    last <- parts[length(parts)]
    if (!grepl("^[0-9]+$", last)) {
        return(NULL)
    }
    cell <- as.integer(last)
    if (is.na(cell)) {
        return(NULL)
    }
    list(
        marker = if (length(parts) == 2L) parts[1L] else NULL,
        cell = cell
    )
}

#' Interactive Ridge Plot From the Ridge Box Inputs
#'
#' @param sce Full \code{SingleCellExperiment}.
#' @param input Shiny input object.
#' @param rs Selected SOM node ids.
#' @param colsUsed Clustering markers (default rows).
#' @param cell Clicked cell index, or \code{NULL}.
#'
#' @return A plotly widget.
#'
#' @keywords internal
.ridgePlotFromInput <- function(sce, input, rs, colsUsed, cell) {
    pick <- input$ridgeClusterSelection
    plotRidgeInteractive(
        sce = sce,
        markers = input$ridgeMarkerSelection %||% colsUsed,
        rs = rs,
        coverage = input$ridgeThreshold %||% 0.99,
        maxCells = max(1L, as.integer(input$ridgeMaxLines %||% 50L)),
        clusterPick = if (!is.null(pick)) as.integer(pick),
        highlightCell = cell,
        clusterMarkers = colsUsed,
        sortRows = isTRUE(input$ridgeSortRows %||% TRUE)
    )
}

#' Node Membership Explanation of the Clicked Cell
#'
#' @param sce Full \code{SingleCellExperiment}.
#' @param somCodes SOM codes matrix.
#' @param colsUsed Clustering markers.
#' @param cell Clicked cell index, or \code{NULL}.
#' @param rs Selected SOM node ids, or \code{NULL}.
#'
#' @return A plotly widget (a hint when no cell is clicked).
#'
#' @keywords internal
.ridgeCellWhyPlot <- function(sce, somCodes, colsUsed, cell, rs) {
    if (is.null(cell) || is.null(rs)) {
        return(plotCellMembership(NULL))
    }
    memb <- .cellMembership(sce, cell, somCodes, colsUsed,
        as.integer(unlist(rs)))
    cd <- SummarizedExperiment::colData(sce)
    lab <- paste("cell", if (is.null(cd$cellId)) cell else cd$cellId[cell])
    plotCellMembership(memb, lab)
}
