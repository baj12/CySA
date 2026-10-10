# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server-reactives-diagnostics.R ----
# Reactives for the Diagnostics, ordered metacluster heatmap, and
# differential-testing panels (see vignettes/catalyst_workflow.Rmd for the
# reference CATALYST/diffcyt figures these reproduce inside CySA).

#' Build the Sample-Diagnostics Reactives
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Called for its side effect of adding reactives to \code{st}.
#'
#' @keywords internal
.buildDiagnosticsReactives <- function(st, input) {
    colorBy <- shiny::reactive({
        cb <- input$diagColorBy
        if (is.null(cb) || identical(cb, "none")) NULL else cb
    })
    diagSce <- shiny::reactive({
        .subsetSceByGroups(st$sce, st$rv$outputList, input$diagGroups)
    })
    st$diagMdsPlot <- shiny::reactive({
        plotSampleMDS(diagSce(), st$colsUsed, color_by = colorBy())
    })
    st$diagExprHeatmapPlot <- shiny::reactive({
        plotSampleExprHeatmap(diagSce(), st$colsUsed)
    })
    st$diagNrsPlot <- shiny::reactive({
        plotNRS(diagSce(), st$colsUsed, color_by = colorBy())
    })
    invisible(st)
}

#' Build the Ordered Cluster Heatmap Reactive
#'
#' @inheritParams .buildDiagnosticsReactives
#'
#' @return Called for its side effect of adding a reactive to \code{st}.
#'
#' @keywords internal
.buildOrderedHeatmapReactives <- function(st, input) {
    st$orderedHeatmapPlot <- shiny::reactive({
        groups <- input$orderedHeatmapGroups %||% names(st$rv$outputList)
        nodeIds <- unique(unlist(
            lapply(st$rv$outputList[groups], as.character),
            use.names = FALSE
        ))
        plotOrderedClusterHeatmap(
            st$metaD[[st$somCodesName]], st$metaD$SOM_stats, st$dend,
            st$dendTable, nodeIds = nodeIds,
            markers = input$orderedHeatmapMarkers %||% st$colsUsed
        )
    })
    invisible(st)
}

#' Initialise the Differential-Testing Result Store
#'
#' The diffcyt GLMM/LMM fit is too slow to run on every input change, so the
#' result is only recomputed when \code{input$runDiffTest} is clicked (see
#' \code{.registerDiffTestObservers()}); this just creates the
#' \code{reactiveVal} holding the latest result.
#'
#' @inheritParams .buildDiagnosticsReactives
#'
#' @return Called for its side effect of adding a reactive value to \code{st}.
#'
#' @keywords internal
.buildDiffTestReactives <- function(st, input) {
    st$diffResult <- shiny::reactiveVal(NULL)
    invisible(st)
}
