# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# outputs_diagnostics.R ----
# Output renderers for the Diagnostics, ordered metacluster heatmap, and
# differential-testing panels (see vignettes/catalyst_workflow.Rmd).

#' Register the Sample-Diagnostics Outputs
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerDiagnosticsOutputs <- function(input, output, st) {
    output$diagMdsPlot <- shiny::renderPlot(st$diagMdsPlot())
    output$diagExprHeatmap <- shiny::renderPlot(st$diagExprHeatmapPlot())
    output$diagNrsPlot <- shiny::renderPlot(st$diagNrsPlot())
    output$diagAmbiguityPlot <- shiny::renderPlot(st$diagAmbiguityPlot())
    invisible(NULL)
}

#' Register the Ordered Cluster Heatmap Output
#'
#' @inheritParams .registerDiagnosticsOutputs
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerOrderedHeatmapOutputs <- function(input, output, st) {
    output$orderedClusterHeatmap <- shiny::renderPlot(st$orderedHeatmapPlot())
    invisible(NULL)
}

#' Register the Differential-Testing Outputs
#'
#' @inheritParams .registerDiagnosticsOutputs
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerDiffTestOutputs <- function(input, output, st) {
    output$diffTestTable <- DT::renderDT(
        options = list(lengthChange = FALSE, scrollX = TRUE), {
            result <- st$diffResult()
            shiny::req(result)
            diffTestTable(result)
        }
    )
    output$diffTestHeatmap <- shiny::renderPlot({
        result <- st$diffResult()
        shiny::req(result)
        tryCatch(
            plotDiffTestHeatmap(result),
            error = function(e) {
                shiny::showNotification(
                    paste("Could not draw the heatmap (the results table",
                        "above is still valid):", conditionMessage(e)),
                    type = "error", duration = NULL
                )
                NULL
            }
        )
    })
    invisible(NULL)
}
