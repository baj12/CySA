# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# outputs-helpers.R ----
# Plot, table and UI builders called by the output renderers of
# outputs_clusterSelector.R. Each takes the server state environment `st`
# and/or plain values, so the renderers themselves stay one-liners.

#' Percentile Axis Limits of a SOM Code Column
#'
#' @param somCodes SOM codes matrix.
#' @param dim Column (marker) name.
#' @param pctl Central fraction to keep, e.g. 0.99.
#'
#' @return Two-element quantile range, or \code{NULL} when \code{dim} is not
#'   a column of \code{somCodes}.
#'
#' @keywords internal
.quantileLimits <- function(somCodes, dim, pctl) {
    if (!dim %in% colnames(somCodes)) {
        return(NULL)
    }
    tailP <- (1 - pctl) / 2
    stats::quantile(somCodes[, dim], probs = c(tailP, 1 - tailP),
        na.rm = TRUE)
}

#' Plotly Cell Scatter Plot with a Selectable Background Grid
#'
#' Wraps the scatter plot of the selected cells with a transparent grid of
#' points so box selections also work in empty areas.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return A plotly widget (empty when no cell is selected).
#'
#' @keywords internal
.scatterPlotly <- function(st, input) {
    rs <- st$rsUsed()
    shiny::req(rs)
    dimSel <- st$dimSelection()
    sampleIds <- input$samples2plot
    shiny::req(dimSel, length(dimSel) >= 1L)
    cd <- SingleCellExperiment::colData(st$sce_subsampled)
    cidIdx <- cd$cluster_id %in% rs & cd$sample_id %in% sampleIds
    if (!any(cidIdx)) {
        return(plotly::plotly_empty())
    }
    pp <- st$scatterPlot()
    if (is.null(pp)) {
        return(plotly::plotly_empty())
    }
    chs <- make.names(dimSel[[1L]]$dims)
    lims <- .computeScatterLimits(
        st$df[cidIdx, chs[1L], drop = TRUE],
        st$df[cidIdx, chs[2L], drop = TRUE],
        input$scatterPercentile %||% 0.99
    )
    gridX <- seq(from = lims$xlim[1L], to = lims$xlim[2L], length.out = 100L)
    gridY <- seq(from = lims$ylim[1L], to = lims$ylim[2L], length.out = 100L)
    quiet_ggplotly(pp, source = "scatterPlot") %>%
        plotly::add_trace(
            x = sort(rep(gridX, 100L)), y = rep(gridY, 100L),
            mode = "markers", type = "scatter",
            fill = "none", fillcolor = "#e763fa", opacity = 0.01
        ) %>%
        plotly::layout(dragmode = "select") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_relayout")
}

#' Plotly SOM 2D Plot Number \code{plotIdxLocal}
#'
#' The last plot follows the active preset pair (\code{activePlot}); with
#' \code{input$showGroups} the nodes are coloured by saved group instead.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#' @param plotIdxLocal Index of the SOM plot output.
#'
#' @return A plotly widget, or \code{NULL}.
#'
#' @keywords internal
.somPlotly <- function(st, input, plotIdxLocal) {
    colorbyGroups <- input$colorbyGroups
    st$selectedUpdate2()
    dimSel <- st$dimSelection()
    rs <- st$rsUsed_d()
    shiny::req(rs)
    st$triggerRedraw()
    plotIdx <- if (plotIdxLocal == st$nPlots) st$activePlot() else plotIdxLocal
    shiny::req(length(dimSel) >= plotIdx)
    dims <- dimSel[[plotIdx]]$dims
    pp1 <- st$getBasePlot(dims[1L], dims[2L], input$somColorVar %||% "n",
        input$somSizeVar %||% "max")
    p3 <- if (isTRUE(input$showGroups)) {
        projectionDf <- st$somProjectionDf()
        shiny::req(projectionDf)
        drawProjection(projectionDf, rs, colorbyGroups = colorbyGroups,
            sce = st$sce, outputList = st$rv$outputList)
    } else {
        ggsomPlot(pp1, plotIdx, rs, dimSelection = dimSel,
            somCodesName = st$somCodesName, sce = st$sce, metaD = st$metaD)
    }
    if (is.null(p3)) {
        return(NULL)
    }
    quiet_ggplotly(p3, source = paste0("somData", plotIdxLocal),
        tooltip = "") %>%
        plotly::layout(showlegend = FALSE, dragmode = "select") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_relayout")
}

#' Main SOM 2D Plot of the Active Axis Pair
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return A plotly widget.
#'
#' @keywords internal
.somMainPlot <- function(st, input) {
    colorbyGroups <- input$colorbyGroups
    st$selectedUpdate2()
    showGroups <- input$showGroups
    nnAgreement <- input$nnAgreement %||% FALSE
    dimSel <- st$dimSelection()
    shiny::req(length(dimSel) >= 1L)
    rs <- st$rsUsed_d()
    shiny::req(rs)
    st$triggerRedraw()
    dims <- dimSel[[1L]]$dims
    pctl <- input$scatterPercentile %||% 0.99
    somCodes <- st$metaD[[st$somCodesName]]
    somPlot(
        st$somBasePlot(), 1L, rs, colorbyGroups, showGroups,
        dimSelection = dimSel, sce = st$sce, metaD = st$metaD,
        outputList = st$rv$outputList,
        projectionDf = if (isTRUE(showGroups) || isTRUE(nnAgreement)) {
            st$somProjectionDf()
        },
        agreement = if (isTRUE(nnAgreement)) st$nnAgreement(),
        xlim = .quantileLimits(somCodes, dims[1L], pctl),
        ylim = .quantileLimits(somCodes, dims[2L], pctl),
        source = "somDataMain"
    )
}

#' Layout of the Static SOM Pair Grid
#'
#' @param current List of marker pairs.
#' @param nCols Number of plots per row (default 3).
#' @param heightPx Height of each plot in pixels (default 240).
#'
#' @return A tag list of fluid rows with one \code{plotOutput} per pair.
#'
#' @keywords internal
.staticSomGridUI <- function(current, nCols, heightPx) {
    nCols <- .asIntegerQuiet(nCols)[1L]
    if (is.na(nCols) || nCols < 1L) nCols <- 3L
    height <- paste0(if (is.null(heightPx)) 240L else heightPx, "px")
    if (length(current) == 0L) {
        return(htmltools::tags$p(
            "No pairs configured. Use the controls above to add pairs.",
            style = "color:#aaa; font-style:italic; padding:12px;"
        ))
    }
    panels <- lapply(seq_along(current), function(i) {
        pair <- current[[i]]
        shiny::column(
            width = 12L %/% nCols,
            htmltools::tags$p(
                paste(pair[1L], "\u00b7", pair[2L]),
                style = paste(
                    "font-size:11px; font-weight:600;",
                    "text-align:center; margin:4px 0 2px; color:#555;"
                )
            ),
            shiny::plotOutput(paste0("staticSomDyn", i), height = height)
        )
    })
    rows <- lapply(seq(1L, length(panels), by = nCols), function(s) {
        do.call(shiny::fluidRow,
            panels[s:min(s + nCols - 1L, length(panels))])
    })
    do.call(htmltools::tagList, rows)
}

#' Static SOM Plot of One Marker Pair
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#' @param pair Character vector of the two markers.
#'
#' @return A ggplot.
#'
#' @keywords internal
.staticSomPlot <- function(st, input, pair) {
    rs <- st$rsUsed_d()
    shiny::req(rs)
    pctl <- input$scatterPercentile %||% 0.99
    pp1 <- st$getBasePlot(pair[1L], pair[2L], input$somColorVar %||% "n",
        input$somSizeVar %||% "max")
    shiny::req(pp1)
    lim1 <- st$channelLimits[[pair[1L]]] %||% c(min = 0, max = 1)
    lim2 <- st$channelLimits[[pair[2L]]] %||% c(min = 0, max = 1)
    dimSel <- list(list(
        dims = pair,
        xlim = c(lim1["min"], lim1["max"]),
        ylim = c(lim2["min"], lim2["max"]),
        xzoom = c(NULL, NULL),
        yzoom = c(NULL, NULL)
    ))
    somCodes <- st$metaD[[st$somCodesName]]
    ggsomPlot(pp1, 1L, rs, dimSel,
        sce = st$sce, metaD = st$metaD,
        xlim = .quantileLimits(somCodes, pair[1L], pctl),
        ylim = .quantileLimits(somCodes, pair[2L], pctl)
    )
}

#' Plotly Version of a Projection Plot
#'
#' @param p3 Projection ggplot (t-SNE, UMAP or PCA).
#' @param source Plotly event source id.
#' @param showLegend Logical; show the legend below the plot.
#'
#' @return A plotly widget with box selection enabled.
#'
#' @keywords internal
.projectionPlotly <- function(p3, source, showLegend) {
    shiny::req(p3)
    ret <- quiet_ggplotly(p3, source = source, tooltip = "text")
    ret <- if (isTRUE(showLegend)) {
        ret %>% plotly::layout(legend = list(
            x = 0, y = -0.2, xanchor = "left", yanchor = "bottom",
            orientation = "h"
        ))
    } else {
        ret %>% plotly::layout(showlegend = FALSE)
    }
    ret %>%
        plotly::layout(dragmode = "select") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_relayout")
}

#' Plotly Dendrogram with Clickable Leaves
#'
#' @param dd Output of \code{.dendPlotlyData()}.
#' @param rs Selected node ids (drawn in red).
#'
#' @return A plotly widget; leaf points carry the node label as
#'   \code{customdata}.
#'
#' @keywords internal
.dendPlotly <- function(dd, rs) {
    labelCols <- ifelse(dd$labels$label %in% as.character(rs), "red",
        "black")
    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = dd$segments,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        ggplot2::geom_point(data = dd$labels, ggplot2::aes(x = x, y = y),
            color = labelCols) +
        ggplot2::geom_text(data = dd$labels,
            ggplot2::aes(x = x, y = y, label = label), vjust = 1, size = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )
    # customdata is injected after ggplotly(); geom_point is trace 2.
    quiet_ggplotly(p, source = "dendPlotly", tooltip = "") %>%
        plotly::style(customdata = dd$labels$label, traces = 2L) %>%
        plotly::layout(dragmode = "select") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_click")
}

#' Format Node Ids in Rows of Forty
#'
#' @param rs Integer node ids.
#'
#' @return A single string, one line per block of 40 node ids.
#'
#' @keywords internal
.formatNodeList <- function(rs) {
    vals <- lapply(seq_len(40), function(x) {
        paste(rs[as.integer(rs / 40) + 1 == x], collapse = ", ")
    })
    paste(unlist(vals), collapse = "\n")
}

#' Numeric Columns of \code{experiment_info}
#'
#' @param metaD \code{metadata(sce)} list with \code{experiment_info}.
#' @param dropCols Column names to exclude.
#'
#' @return A data frame of the numeric columns, sample ids as row names.
#'
#' @keywords internal
.numericExperimentInfo <- function(metaD, dropCols = character(0)) {
    info <- metaD$experiment_info
    numCols <- unlist(lapply(info, is.numeric), use.names = FALSE)
    expInfo <- info[, numCols, drop = FALSE]
    expInfo <- expInfo[, !colnames(expInfo) %in% dropCols, drop = FALSE]
    eI <- apply(expInfo, 2, as.numeric) %>% as.data.frame()
    rownames(eI) <- info$sample_id
    eI
}
