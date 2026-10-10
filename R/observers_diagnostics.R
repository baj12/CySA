# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# observers_diagnostics.R ----
# Observer for the diffcyt differential-testing panel (see
# vignettes/catalyst_workflow.Rmd).

#' Register the Differential-Testing Observer
#'
#' Runs a diffcyt DA/DS test on the current saved groups
#' (\code{st$rv$outputList}) when "Run test" is clicked, and stores the
#' result in \code{st$diffResult}. The fit is gated behind an action button
#' rather than reactive, since GLMM/LMM fitting is too slow to rerun on
#' every input change.
#'
#' @param input,output,session Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerDiffTestObservers <- function(input, output, session, st) {
    shiny::observeEvent(input$runDiffTest, {
        outputList <- st$rv$outputList[intersect(
            input$diffGroups, names(st$rv$outputList)
        )]
        shiny::req(
            length(outputList) >= 2, input$diffFixedCol,
            nzchar(input$diffFixedCol)
        )
        if (identical(input$diffAnalysisType, "DA") &&
            identical(input$diffRandomCol, "none")) {
            shiny::showNotification(paste(
                "DA (diffcyt-DA-GLMM) requires a random effect, e.g.",
                "sample_id -- lme4::glmer() cannot fit a model with no",
                "random-effect term."
            ), type = "error", duration = NULL)
            return(invisible(NULL))
        }
        if (identical(input$diffAnalysisType, "DS") &&
            identical(input$diffRandomCol, "sample_id")) {
            shiny::showNotification(paste(
                "DS test: sample_id as the random effect is confounded",
                "with the residual (each cluster-sample combination has",
                "one median value); pick 'none' or another column (e.g.",
                "patient_id)."
            ), type = "error", duration = NULL)
            return(invisible(NULL))
        }
        markers <- if (identical(input$diffAnalysisType, "DS")) {
            shiny::req(length(input$diffMarkers) > 0)
            input$diffMarkers
        } else {
            NULL
        }
        result <- withCallingHandlers(
            tryCatch(
                runDiffTest(
                    sce = st$sce, outputList = outputList,
                    experimentInfo = st$metaD$experiment_info,
                    fixedCol = input$diffFixedCol,
                    randomCol = input$diffRandomCol,
                    analysisType = input$diffAnalysisType,
                    markers = markers
                ),
                error = function(e) {
                    shiny::showNotification(
                        paste("Differential test failed:",
                            conditionMessage(e)),
                        type = "error", duration = NULL
                    )
                    NULL
                }
            ),
            warning = function(w) {
                shiny::showNotification(conditionMessage(w),
                    type = "warning", duration = 10)
                invokeRestart("muffleWarning")
            }
        )
        if (!is.null(result)) st$diffResult(result)
    })
    invisible(NULL)
}
