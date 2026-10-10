# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# interface_diagnostics.R ----
# UI builders for the Diagnostics, ordered metacluster heatmap, and
# differential-testing panels, added on top of the CATALYST/diffcyt
# reference workflow reproduced in vignettes/catalyst_workflow.Rmd.

#' Column-Name Choices From experiment_info
#'
#' @param info \code{experiment_info} data frame.
#'
#' @return A list with \code{factor} and \code{numeric} column-name vectors,
#'   each prefixed with \code{"none"}.
#'
#' @keywords internal
.experimentInfoChoices <- function(info) {
    factCols <- colnames(info)[unlist(lapply(info, is.factor),
        use.names = FALSE)]
    numCols <- colnames(info)[unlist(lapply(info, is.numeric),
        use.names = FALSE)]
    list(factor = c("none", factCols), numeric = c("none", numCols))
}

#' Build Diagnostics UI
#'
#' Sample-level QC panel: pseudobulk MDS, a per-sample marker-expression
#' heatmap, and marker non-redundancy scores (NRS) -- CySA's equivalents of
#' \code{CATALYST::pbMDS()}/\code{plotExprHeatmap()}/\code{plotNRS()} (see
#' \code{vignettes/catalyst_workflow.Rmd}).
#'
#' @param metaD \code{metadata(sce)} list.
#' @param outputList Named list of cluster groupings, used to restrict the
#'   plots to a selection of clusters/metaclusters.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildDiagnosticsBox <- function(metaD, outputList) {
    choices <- .experimentInfoChoices(metaD$experiment_info)
    colorChoices <- union(choices$factor, choices$numeric)
    .primaryBox(
        "Sample diagnostics", id = "diagnosticsBox", infoId = "diagnostics",
        shiny::fluidRow(
            .tipColumn(4, "diagGroups", shiny::selectInput("diagGroups",
                "Clusters/groups to include", choices = names(outputList),
                selected = names(outputList), multiple = TRUE,
                selectize = TRUE)),
            .tipColumn(4, "diagColorBy", shiny::selectInput("diagColorBy",
                "Colour samples by", choices = colorChoices,
                selected = "none"))
        ),
        shiny::tabsetPanel(
            shiny::tabPanel("Sample MDS",
                .resizable(shiny::plotOutput("diagMdsPlot"))),
            shiny::tabPanel("Expression heatmap",
                .resizable(shiny::plotOutput("diagExprHeatmap"))),
            shiny::tabPanel("Marker NRS",
                .resizable(shiny::plotOutput("diagNrsPlot")))
        )
    )
}

#' Build Ordered Cluster Heatmap UI
#'
#' Dendrogram-ordered marker-by-cluster heatmap with a cell-count bar --
#' CySA's equivalent of \code{CATALYST::plotMultiHeatmap(bars = TRUE)},
#' reusing the existing SOM dendrogram (see
#' \code{vignettes/catalyst_workflow.Rmd}).
#'
#' @param outputList Named list of cluster groupings, used to restrict the
#'   heatmap to a selection of clusters/metaclusters.
#' @param markers Marker choices offered for the heatmap's columns.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildOrderedHeatmapBox <- function(outputList, markers) {
    .primaryBox(
        "Ordered cluster heatmap", id = "orderedHeatmapBox",
        infoId = "orderedHeatmap",
        shiny::fluidRow(
            .tipColumn(4, "orderedHeatmapGroups", shiny::selectInput(
                "orderedHeatmapGroups", "Clusters/groups to show",
                choices = names(outputList), selected = names(outputList),
                multiple = TRUE, selectize = TRUE)),
            .tipColumn(4, "orderedHeatmapMarkers", shiny::selectInput(
                "orderedHeatmapMarkers", "Markers to show",
                choices = markers, selected = markers, multiple = TRUE,
                selectize = TRUE))
        ),
        shiny::fluidRow(shiny::column(width = 12,
            .resizable(shiny::plotOutput("orderedClusterHeatmap",
                height = "700px"))))
    )
}

#' Build Differential Testing UI
#'
#' diffcyt GLMM/LMM differential abundance (DA) and differential state (DS)
#' tests on CySA's own saved cluster groups, rather than on raw SOM nodes
#' (see \code{vignettes/catalyst_workflow.Rmd} for the reference workflow on
#' CATALYST metaclusters).
#'
#' @param metaD \code{metadata(sce)} list.
#' @param outputList Named list of cluster groupings; the test compares a
#'   user-chosen selection of these clusters/metaclusters.
#' @param markers Marker choices offered for the DS test.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildDiffTestBox <- function(metaD, outputList, markers) {
    choices <- .experimentInfoChoices(metaD$experiment_info)
    factChoices <- setdiff(choices$factor, "none")
    randomChoices <- union("sample_id", factChoices)
    # "selected" mirrors whatever nodes were just picked/named and almost
    # always overlaps a named group, which a differential test cannot
    # tolerate (see .checkGroupsDisjoint()); default it off.
    defaultGroups <- setdiff(names(outputList), "selected")
    .primaryBox(
        "Differential testing (diffcyt)", id = "diffTestBox",
        infoId = "diffTest",
        htmltools::tags$p(
            "Tests the chosen saved groups for differential abundance (DA)",
            "or differential marker state (DS) across samples, using",
            "diffcyt's GLMM/LMM methods. The compared groups must not",
            "overlap in SOM nodes -- each node must belong to exactly one",
            "of them."
        ),
        shiny::fluidRow(
            .tipColumn(6, "diffGroups", shiny::selectInput("diffGroups",
                "Clusters/groups to compare (>= 2, non-overlapping)",
                choices = names(outputList), selected = defaultGroups,
                multiple = TRUE, selectize = TRUE)),
            .tipColumn(6, "diffMarkers", shiny::selectInput("diffMarkers",
                "Markers to test (DS only)", choices = markers,
                selected = markers[1], multiple = TRUE, selectize = TRUE))
        ),
        shiny::fluidRow(
            .tipColumn(4, "diffAnalysisType", shiny::radioButtons(
                "diffAnalysisType", "Analysis type",
                choices = c("Abundance (DA)" = "DA", "State (DS)" = "DS"),
                selected = "DA")),
            .tipColumn(4, "diffFixedCol", shiny::selectInput("diffFixedCol",
                "Fixed effect (2-level factor)", choices = factChoices)),
            .tipColumn(4, "diffRandomCol", shiny::selectInput("diffRandomCol",
                "Random effect", choices = c("none", randomChoices),
                selected = "sample_id"))
        ),
        htmltools::tags$p(shiny::tags$em(
            "DA (diffcyt-DA-GLMM) requires a random effect -- normally",
            "sample_id. For DS, do not use sample_id as the random effect:",
            "each cluster-sample combination already contributes one",
            "median value, so a sample_id random intercept is perfectly",
            "confounded with the residual and the fit silently returns NA",
            "p-values for every row; pick 'none' or a column such as",
            "patient_id instead."
        )),
        shiny::fluidRow(shiny::column(width = 12,
            .withTip("runDiffTest", shiny::actionButton("runDiffTest",
                "Run test")))),
        shiny::fluidRow(shiny::column(width = 12,
            .resizable(DT::DTOutput("diffTestTable")))),
        shiny::fluidRow(shiny::column(width = 12,
            .resizable(shiny::plotOutput("diffTestHeatmap",
                height = "600px"))))
    )
}
