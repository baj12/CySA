# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server-reactive-helpers.R ----
# Reactive helper functions extracted from the clusterSelector() server closure.
# These functions build the dimension-reduction reactives, selection observers,
# and zoom observers used by the app. They are kept in a separate file so the
# main server reads as pure wiring.

#' Build Dimension Reduction Reactives
#'
#' Creates the t-SNE, UMAP, and PCA reactives used by the clusterSelector app.
#'
#' @param input Shiny input object.
#' @param metaD \code{metadata(sce)} list.
#' @param sce Full \code{SingleCellExperiment}.
#' @param somCodesName Name of the SOM codes metadata slot.
#'
#' @return A named list of three reactives: \code{tsne}, \code{umap},
#'   \code{pca}.
#'
#' @keywords internal
.buildDimRedReactives <- function(input, metaD, sce, somCodesName) {
    tsne <- shiny::reactive({
        tsneFunc(input$dimRedSelection, input$perplexity, sce, somCodesName)
    }) %>%
        shiny::bindCache(input$dimRedSelection, input$perplexity) %>%
        shiny::debounce(1000)

    umap <- shiny::reactive({
        cols <- input$dimRedSelection
        pumap <- umap::umap.defaults
        pumap$n_neighbors <- input$n_neighbors
        umap::umap(metaD[[somCodesName]][, cols, drop = FALSE], config = pumap)
    }) %>%
        shiny::bindCache(input$dimRedSelection, input$n_neighbors) %>%
        shiny::debounce(1000)

    pca <- shiny::reactive({
        cols <- input$dimRedSelection
        stats::prcomp(t(metaD[[somCodesName]][, cols, drop = FALSE]),
            scale = FALSE, rank. = 2)
    }) %>%
        shiny::bindCache(input$dimRedSelection) %>%
        shiny::debounce(1000)

    list(tsne = tsne, umap = umap, pca = pca)
}


#' Build Selection Observer
#'
#' Registers a \code{plotly_selected} observer for a single source. When a
#' selection event fires, it updates \code{rsUsed} according to the current
#' selection mode.
#'
#' @param sourceId Plotly source id.
#' @param input Shiny input object.
#' @param rsUsed Reactive value holding the currently selected nodes.
#' @param rsUsed_d Debounced reactive read as the selection baseline before
#'   applying \code{inputSelect}.
#' @param inputSelect Function that maps event data to new selected ids.
#' @param verbose Logical indicating whether to print debug messages.
#'
#' @keywords internal
.buildSelectionObserver <- function(sourceId, input, rsUsed, rsUsed_d,
    inputSelect, verbose) {
    # Callers register observers in a loop: evaluate the source id now, not
    # when the observer first fires.
    force(sourceId)
    shiny::observeEvent(
        .safeEventData("plotly_selected", source = sourceId,
            verbose = verbose),
        {
            d <- .safeEventData("plotly_selected", source = sourceId,
                verbose = verbose)
            if (is.null(d)) {
                return(NULL)
            }
            rs <- shiny::isolate(rsUsed_d())
            if (is.null(rs)) {
                return(NULL)
            }
            d <- inputSelect(d, rs, shiny::isolate(input$selectMode))
            shiny::isolate(rsUsed(d))
        }
    )
}

#' Build SOM Data Selection Observers
#'
#' Loops over the 2D SOM plot outputs and registers selection observers.
#'
#' @param nPlots Number of 2D SOM plots.
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param rsUsed Reactive value holding the currently selected nodes.
#' @param rsUsed_d Debounced reactive read as the selection baseline.
#' @param activePlot Reactive value tracking the active plot index.
#' @param inputSelect Function that maps event data to new selected ids.
#' @param verbose Logical indicating whether to print debug messages.
#'
#' @keywords internal
.buildSOMDataObservers <- function(nPlots, input, output, rsUsed, rsUsed_d,
    activePlot, inputSelect, verbose) {
    lapply(seq_len(nPlots), function(i) {
        src <- paste0("somData", i)
        shiny::observeEvent(
            .safeEventData("plotly_selected", source = src,
                verbose = verbose),
            {
                if (isTRUE(verbose)) message("som", i, " touched")
                activePlot(i)
                rs <- shiny::isolate(rsUsed_d())
                shiny::req(rs)
                d <- .safeEventData("plotly_selected", source = src,
                    verbose = verbose)
                if (is.null(d)) {
                    return(NULL)
                }
                d <- inputSelect(d, rs, shiny::isolate(input$selectMode))
                shiny::isolate(rsUsed(d))
            }
        )
    })
}

#' Build Zoom Observers for SOM Data Plots
#'
#' Registers one \code{plotly_relayout} observer per 2D SOM plot.
#'
#' @param nPlots Number of 2D SOM plots.
#' @param input Shiny input object.
#' @param zoomFunc Function that applies zoom state to \code{dimSelection}.
#' @param verbose Logical indicating whether to print debug messages.
#'
#' @keywords internal
.buildZoomObservers <- function(nPlots, input, zoomFunc, verbose) {
    lapply(seq_len(nPlots), function(plotIdx) {
        shiny::observe({
            zoom <- .safeEventData(
                verbose = verbose, "plotly_relayout",
                source = paste0("somData", plotIdx)
            )
            if (is.null(zoom)) {
                return(NULL)
            }
            zoomFunc(zoom, plotIdx)
        })
    })
}
