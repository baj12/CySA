# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# outputs-helpers-stats.R ----
# Statistics tables, raster / FlowSOM output builders and the PDF plot
# registry used by the output renderers of outputs_clusterSelector.R.

#' Per-Sample Cell Counts of the Selection
#'
#' One row with the per-sample cell counts of the selected nodes, optionally
#' a comparison row (a numeric \code{experiment_info} column or a saved
#' group, \code{input$compareStatsTo}) and the counts of one single node
#' (\code{input$singleNode}).
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return A \code{data.table}.
#'
#' @keywords internal
.cellCountsTable <- function(st, input) {
    cst <- input$compareStatsTo
    sN <- as.integer(input$singleNode)
    rs <- st$rsUsed()
    shiny::req(rs)
    outputList <- st$rv$outputList
    cpt <- st$clusterPatientTable
    shiny::req(cpt)
    nodeCounts <- function(nodes) {
        nodes <- .validNodeCols(nodes, cpt)
        data.table::as.data.table(t(rowSums(cpt[, nodes, drop = FALSE])))
    }
    rSums <- nodeCounts(rs)
    colnames(rSums) <- rownames(cpt)
    if (!identical(cst, "none")) {
        ctStats <- NULL
        if (cst %in% colnames(st$metaD$experiment_info)) {
            eI <- .numericExperimentInfo(st$metaD,
                c("sample_nr", "sample_id", "sample"))
            if (any(is.na(eI))) {
                stop("NAs produced when converting experiment_info to numeric")
            }
            ctStats <- data.table::as.data.table(
                t(eI[colnames(rSums), cst, drop = FALSE])
            )
        } else if (cst %in% names(outputList)) {
            ctStats <- nodeCounts(outputList[[cst]])
        }
        if (!is.null(ctStats)) {
            colnames(ctStats) <- colnames(rSums)
            rSums <- data.table::rbindlist(list(rSums, ctStats))
            rownames(rSums) <- c("selection", cst)
        }
    }
    rn <- rownames(rSums)
    ctStats <- nodeCounts(sN)
    colnames(ctStats) <- rownames(cpt)
    rSums <- data.table::rbindlist(list(rSums, ctStats))
    rownames(rSums) <- c(rn, paste("SOM node", sN))
    rSums
}

#' Per-Sample Percentages of the Selection
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Formatted percentages per sample (\code{noquote}).
#'
#' @keywords internal
.cellPercentages <- function(st, input) {
    outputList <- st$rv$outputList
    rs <- st$rsUsed()
    shiny::req(rs, st$clusterPatientTable)
    # .numericExperimentInfo() re-attaches sample_id row names after
    # apply()/as.numeric(); building eI inline here previously lost them,
    # so indexing by sample id below threw "subscript out of bounds".
    eI <- .numericExperimentInfo(st$metaD)
    rSums <- .computeRelativeCounts(
        st$clusterPatientTable, rs, input$relativeTo, eI, outputList
    )
    names(rSums) <- rownames(st$clusterPatientTable)
    noquote(formatC(signif(rSums, digits = 2), digits = 2, format = "fg",
        flag = "#"))
}

#' Two-Sample t-Test of the Selection Between Sample Groups
#'
#' Compares the per-sample cell counts of the selected nodes (optionally
#' relative to \code{input$relativeTo}) between the samples of
#' \code{group1} and \code{group2}.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return An \code{htest} object, or a message string when there is not
#'   enough data.
#'
#' @keywords internal
.ttestResult <- function(st, input) {
    groupsVar <- input$groupsVar
    empty <- is.null(groupsVar) || length(groupsVar) == 0 ||
        !groupsVar %in% colnames(st$metaD$experiment_info)
    grpInp <- st$groupsInput()
    empty <- empty || purrr::is_empty(grpInp) ||
        any(unlist(lapply(grpInp, purrr::is_empty)))
    rs <- st$rsUsed()
    shiny::req(rs)
    if (empty || length(rs) < 1) {
        return("no data")
    }
    cpt <- st$clusterPatientTable
    rSums <- if (input$relativeTo == "none") {
        rep(1, nrow(cpt))
    } else {
        .computeRelativeCounts(cpt, rs, input$relativeTo,
            .numericExperimentInfo(st$metaD), st$rv$outputList
        )[rownames(cpt)] / 100
    }
    names(rSums) <- rownames(cpt)
    cD <- SingleCellExperiment::colData(st$sce)
    cellCounts <- cD[cD$cluster_id %in% rs, ] %>%
        dplyr::as_tibble() %>%
        dplyr::group_by(.data$sample_id) %>%
        dplyr::count()
    x <- .groupRates(cellCounts, grpInp$group1, rSums)
    y <- .groupRates(cellCounts, grpInp$group2, rSums)
    shiny::req(x, y)
    if (stats::sd(x) == 0 && stats::sd(y) == 0) {
        return("not enough data")
    }
    stats::t.test(x, y)
}

#' Per-Sample Rates of One Sample Group
#'
#' @param cellCounts Tibble with \code{sample_id} and count \code{n}.
#' @param samples Sample ids of the group.
#' @param rSums Named per-sample denominators.
#'
#' @return Numeric vector \code{n / rSums} for the samples of the group.
#'
#' @keywords internal
.groupRates <- function(cellCounts, samples, rSums) {
    x <- cellCounts %>%
        dplyr::filter(.data$sample_id %in% samples) %>%
        dplyr::ungroup()
    x$n / rSums[x$sample_id]
}

#' SOM Raster Heatmap With the Selected Nodes Marked
#'
#' @param baseRasterGgplot Pre-built raster ggplot.
#' @param xy Data frame with the \code{x}, \code{y} grid positions of the
#'   selected nodes.
#'
#' @return A ggplot.
#'
#' @keywords internal
.rasterWithSelection <- function(baseRasterGgplot, xy) {
    baseRasterGgplot +
        ggplot2::geom_point(
            data = xy, ggplot2::aes(x = x, y = y),
            color = "red", size = 1, inherit.aes = FALSE
        ) +
        ggplot2::geom_point(
            data = xy, ggplot2::aes(x = x, y = y),
            color = "red", shape = 3, size = 1, inherit.aes = FALSE
        )
}

#' Sized Plot Output for the FlowSOM Marker Pies
#'
#' @param n Number of selected nodes.
#' @param ncol Pies per row.
#' @param size Pixel size of one pie.
#'
#' @return A \code{plotOutput}, or \code{NULL}.
#'
#' @keywords internal
.flowSOMPieUI <- function(n, ncol, size) {
    if (is.null(ncol) || is.null(size) || n < 1L) {
        return(NULL)
    }
    # Limit the dimensions to keep the PNG well below its 50000 px limit.
    ncol <- min(ncol, 20L)
    size <- min(size, 300L)
    maxDim <- 8000L
    shiny::plotOutput(
        "flowSOMPie",
        width = paste0(min(ncol * size, maxDim), "px"),
        height = paste0(min(ceiling(n / ncol) * size, maxDim), "px")
    )
}

#' Plotly FlowSOM Star Plot With Selectable Nodes
#'
#' @param fsom A \code{FlowSOM} object.
#'
#' @return A plotly widget; invisible node points carry the node id as
#'   \code{customdata}.
#'
#' @keywords internal
.flowSOMStarsPlotly <- function(fsom) {
    shiny::req(fsom)
    ret <- FlowSOM::PlotStars(fsom, backgroundValues = fsom$metaclustering,
        list_insteadof_ggarrange = TRUE)
    coords <- fsom$MST$l
    coordDf <- data.frame(x = coords[, 1], y = coords[, 2],
        id = seq_len(nrow(coords)))
    pClick <- ret$tree +
        ggplot2::geom_point(data = coordDf, ggplot2::aes(x = x, y = y),
            alpha = 0, size = 8, inherit.aes = FALSE)
    # The invisible node points are the last trace.
    lastTrace <- length(quiet_ggplotly(pClick)$x$data)
    quiet_ggplotly(pClick, source = "flowSOMStars", tooltip = "") %>%
        plotly::style(customdata = coordDf$id, traces = lastTrace) %>%
        plotly::layout(dragmode = "select") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_click")
}

#' Cluster Picker of the Marker Ridge Box
#'
#' Offers the selected SOM nodes that contain cells. Individual density
#' lines are meant for a handful of nodes, so the picker pre-selects the
#' five largest nodes and accepts at most \code{.ridgeMaxLineClusters}.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id}.
#' @param rs Selected node ids.
#' @param current Currently picked nodes, or \code{NULL}.
#'
#' @return A tag list with the picker, or \code{NULL} with five or fewer
#'   selected nodes.
#'
#' @keywords internal
.ridgeClusterUI <- function(sce, rs, current) {
    cid <- as.integer(as.character(
        SingleCellExperiment::colData(sce)$cluster_id
    ))
    sizes <- table(cid[cid %in% as.integer(rs)])
    if (length(sizes) <= 5L) {
        return(NULL)
    }
    clusterIds <- sort(as.integer(names(sizes)))
    largest <- sort(as.integer(names(sort(-sizes)[seq_len(5L)])))
    shiny::tagList(shiny::div(
        title = .clusterSelectorTooltips[["ridgeClusterSelection"]],
        shiny::selectizeInput(
            "ridgeClusterSelection", "Clusters with individual lines",
            choices = as.character(clusterIds),
            selected = current %||% as.character(largest),
            multiple = TRUE,
            options = list(maxItems = .ridgeMaxLineClusters)
        )
    ))
}

#' Maximum Number of Individual Node Lines in the Ridge Plots
#'
#' @keywords internal
.ridgeMaxLineClusters <- 20L

#' Plots Written to the PDF Download
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#' @param rs Selected node ids.
#'
#' @return Named list of functions (or \code{NULL}) returning plots.
#'
#' @keywords internal
.pdfPlotRegistry <- function(st, input, rs) {
    sce <- st$sce
    upsetSel <- input$upsetSelection
    violinSel <- input$violinSelection %||% st$colsUsed
    outList <- st$rv$outputList
    if (length(upsetSel) < 3) upsetSel <- names(outList)
    list(
        countBar = st$countBarPlot, percentBar = st$PercentBarPlot,
        tsne = st$tsnePlot, umap = st$umapPlot, pca = st$pcaPlot,
        scatter = st$scatterPlot, flowSOMPie = st$flowSOMPiePlot,
        flowSOMStars = if (!is.null(st$fsom)) {
            function() {
                FlowSOM::PlotStars(st$fsom,
                    backgroundValues = st$fsom$metaclustering)$tree
            }
        },
        violin1 = function() {
            plotViolinFunc(sce, st$somCodesName, upsetSel, outList,
                violinSel)
        },
        violin2 = function() {
            plotViolin2Func(sce, st$somCodesName, violinSel, upsetSel,
                outList)
        },
        ridgeTrace = function() .pdfRidgeTrace(st, input, rs),
        upset = function() upsetPlotFunc(upsetSel, outList, sce)
    )
}

#' Static Marker Ridge Plot for the PDF Download
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#' @param rs Selected node ids.
#'
#' @return A ggplot.
#'
#' @keywords internal
.pdfRidgeTrace <- function(st, input, rs) {
    ridgeClusters <- input$ridgeClusterSelection
    plotRidgeTraceFunc(
        sce = st$sce,
        somCodesName = st$somCodesName,
        markers = input$ridgeMarkerSelection %||% st$colsUsed,
        rs = as.integer(unlist(rs)),
        # The input is a band coverage (fraction of cells inside); the
        # static plot takes an sd multiple, so use the normal quantile with
        # the same two-sided coverage (0.99 -> 2.58 sd).
        threshold = stats::qnorm(
            (1 + min(input$ridgeThreshold %||% 0.99, 0.999999)) / 2
        ),
        maxLines = max(1L, as.integer(input$ridgeMaxLines %||% 50L)),
        maxClusters = 5L,
        clusterPick = if (!is.null(ridgeClusters)) {
            as.integer(ridgeClusters)
        }
    )
}
