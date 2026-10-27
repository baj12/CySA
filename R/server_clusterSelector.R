# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server_clusterSelector.R ----
# Server factory for the clusterSelector Shiny app.
# This file builds the `server` function that is passed to shiny::shinyApp().
#
# All per-session state lives in one environment `st`: static data derived
# from the inputs, reactive values, and reactives. Builders add to `st` in
# stages; reactive bodies read `st$<name>` lazily when they run, so the
# order in which builders define them does not matter. At the end the
# bindings are copied into the server function's own environment, so
# shiny::testServer() expressions (and debugging) see `rsUsed`, `rv`,
# `dimSelection`, ... by name.

#' Build clusterSelector Server
#'
#' Creates the Shiny server function for the clusterSelector app.
#'
#' @param sce Full \code{SingleCellExperiment}.
#' @param sce_subsampled Subsampled \code{SingleCellExperiment}.
#' @param outputList Named list of cluster groupings.
#' @param colTree Optional collapsible tree object.
#' @param dList List of marker pairs for 2D plots.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param nPlots Number of 2D SOM plots to display.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somRasterObj Raster object for SOM visualization.
#' @param fsom Optional \code{FlowSOM} object.
#' @param env Environment used to store mutable state (legacy argument).
#' @param verbose Logical indicating whether to show detailed runtime messages.
#'
#' @importFrom tidyselect all_of
#'
#' @return A \code{shiny} server function.
#'
#' @keywords internal
.buildClusterSelectorServer <- function(
        sce, sce_subsampled, outputList, colTree, dList, dend, dendTable,
        clusterPatientTable, somCodesName, nPlots, somRasterData,
        somRasterObj, fsom, env, verbose) {
    args <- as.list(environment())
    function(input, output, session) {
        # nocov start
        st <- .initServerState(args, input, session)
        .buildDimRedPlotReactives(st, input)
        .buildSelectionPlotReactives(st, input)
        .buildStatsReactives(st, input)
        .buildGroupPlotReactives(st, input)
        .buildDiagnosticsReactives(st, input)
        .buildNodePropsReactives(st, input)
        .buildOrderedHeatmapReactives(st, input)
        .buildDiffTestReactives(st, input)
        st$ridgeServer <- .buildRidgeServer(
            input, output, st$sce,
            somCodes = st$metaD[[st$somCodesName]], colsUsed = st$colsUsed,
            rsUsed = st$rsUsed, selectedUpdate = st$selectedUpdate,
            refresh = function() {
                input$applyName
                input$rmGrp
            },
            verbose = st$verbose
        )
        .buildVerboseLoggers(st, input)
        # Expose the session state by name in this environment.
        list2env(as.list(st, all.names = TRUE), environment())
        .registerClusterSelectorObservers(input, output, session, st)
        .registerClusterSelectorOutputs(input, output, session, st)
        # nocov end
    }
}

#' Initialise the Per-Session Server State
#'
#' Creates the state environment holding the static inputs, derived static
#' data (row/column names, subsampled data frame, per-channel limits, the
#' raster background) and the reactive values of one app session.
#'
#' @param args Named list of the arguments of
#'   \code{.buildClusterSelectorServer()}.
#' @param input,session Shiny objects.
#'
#' @return The state environment.
#'
#' @keywords internal
.initServerState <- function(args, input, session) {
    st <- list2env(args, envir = new.env(parent = emptyenv()))
    st$rnSCE <- rownames(st$sce)
    st$metaD <- S4Vectors::metadata(st$sce)
    st$colsUsed <- st$metaD$map$colsUsed
    ## Reactive values
    st$rsUsed <- shiny::reactiveVal(c(1))
    st$triggerRedraw <- shiny::reactiveVal(1)
    st$selectedPoints <- shiny::reactiveVal(NULL)
    st$activePlot <- shiny::reactiveVal(1)
    st$dListRV <- shiny::reactiveVal(st$dList)
    st$selectedUpdate2 <- shiny::reactiveVal(value = 0)
    st$selectedUpdate <- shiny::reactiveVal(value = 0)
    st$dimSelection <- shiny::reactiveVal(list())
    st$rv <- shiny::reactiveValues(outputList = st$outputList)
    st$rsUsed_d <- st$rsUsed %>% shiny::debounce(1500)
    st$sample2PlotDb <- shiny::reactive(input$samples2plot) %>%
        shiny::debounce(500)
    ## Static / cached data
    st$df <- as.data.frame(t(
        SummarizedExperiment::assays(st$sce_subsampled)[[1]]
    ))
    colnames(st$df) <- make.names(rownames(st$sce_subsampled))
    st$dfPlot <- shiny::reactiveVal(st$df)
    st$sce_subsampledRN <- rownames(st$sce_subsampled)
    st$sceRN <- rownames(st$sce)
    st$sceCN <- colnames(st$sce)
    st$channelLimits <- .computeChannelLimits(st$sce)
    st$baseRasterGgplot <- .buildBaseRasterGgplot(
        st$somRasterData, st$colsUsed
    )
    st$getBasePlot <- .makeBasePlotCache(st)
    st$updatedoutputList <- function() {
        .updateOutputListInputs(session, input, st$rv$outputList, st$metaD)
    }
    st$zoomFunc <- function(zoom, plotIdx) .applyZoom(st, zoom, plotIdx)
    st
}

#' Per-Channel Axis Limits
#'
#' @param sce A \code{SingleCellExperiment}; the first assay is used.
#'
#' @return Named list of \code{c(min = , max = )} per row of \code{sce}.
#'
#' @keywords internal
.computeChannelLimits <- function(sce) {
    assayMat <- SummarizedExperiment::assays(sce)[[1]]
    lims <- lapply(rownames(sce), function(ch) {
        vals <- assayMat[ch, ]
        c(min = min(vals), max = max(vals))
    })
    names(lims) <- rownames(sce)
    lims
}

#' Cached SOM Base Plot Factory
#'
#' Returns a function building the SOM scatter base plot for a channel pair,
#' colour and size variable. Results are cached per session so that only
#' selection changes do not re-run \code{plotSOMScatter()}.
#'
#' @param st Server state environment.
#'
#' @return A function \code{(ch1, ch2, colorVar, sizeVar)}.
#'
#' @keywords internal
.makeBasePlotCache <- function(st) {
    cache <- new.env(parent = emptyenv())
    function(ch1, ch2, colorVar, sizeVar) {
        key <- paste(ch1, ch2, colorVar, sizeVar, sep = "\x01")
        if (!exists(key, envir = cache, inherits = FALSE)) {
            cache[[key]] <- plotSOMScatter(
                x = st$sce, chs = c(ch1, ch2), pointSize = sizeVar,
                color_by = colorVar, zeros = TRUE,
                xRN = st$sceRN, xCN = st$sceCN
            )
        }
        cache[[key]]
    }
}

#' Apply a Plotly Zoom Event to the Dimension Selection
#'
#' Stores the zoomed axis ranges of plot \code{plotIdx} in
#' \code{st$dimSelection}; zooming out past the stored range resets it.
#'
#' @param st Server state environment.
#' @param zoom \code{plotly_relayout} event data, or \code{NULL}.
#' @param plotIdx Index of the SOM plot.
#'
#' @return Called for its side effect.
#'
#' @keywords internal
.applyZoom <- function(st, zoom, plotIdx) {
    dimSel <- st$dimSelection()
    shiny::req(dimSel)
    if (!all(!is.null(zoom), "xaxis.range[0]" %in% names(zoom),
        na.rm = TRUE)) {
        return(invisible(NULL))
    }
    rezoom <- all(c(
        dimSel[[plotIdx]]$xzoom[1] > zoom$`xaxis.range[0]`,
        !is.null(dimSel[[plotIdx]]$xzoom[1])
    ), na.rm = TRUE)
    if (rezoom) {
        dimSel[[plotIdx]]$xzoom <- c(NULL, NULL)
        dimSel[[plotIdx]]$yzoom <- c(NULL, NULL)
    } else {
        dimSel[[plotIdx]]$xzoom <-
            c(zoom$`xaxis.range[0]`, zoom$`xaxis.range[1]`)
        dimSel[[plotIdx]]$yzoom <-
            c(zoom$`yaxis.range[0]`, zoom$`yaxis.range[1]`)
    }
    st$dimSelection(dimSel)
    shiny::isolate(st$triggerRedraw(st$triggerRedraw() + 1))
}

#' Build the Verbose Logging Observers
#'
#' Only active when \code{st$verbose} is \code{TRUE}: keeps the
#' dimension-reduction and SOM base plot reactives warm so their computation
#' shows up in the log.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Called for its side effect.
#'
#' @keywords internal
.buildVerboseLoggers <- function(st, input) {
    if (!isTRUE(st$verbose)) {
        return(invisible(NULL))
    }
    shiny::observe(st$tsne()) %>%
        shiny::bindCache(st$dimRedSelection(), input$perplexity)
    shiny::observe(st$umap()) %>%
        shiny::bindCache(input$dimRedSelection, input$n_neighbors)
    shiny::observe(st$pca()) %>%
        shiny::bindCache(input$dimRedSelection)
    shiny::observe(st$somBasePlot()) %>%
        shiny::bindCache(
            st$activeDims(), input$somColorVar, input$somSizeVar
        )
    invisible(NULL)
}
