# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# interface_clusterSelector.R ----
# UI construction helpers for clusterSelector(). Splitting the dashboard into
# small builder functions makes the main factory easier to read and test.


# Central tooltip texts -------------------------------------------------------
# Named character vector used to keep hover-over help consistent across all UI
# builder functions.  Names correspond to Shiny inputIds so that adding a new
# control only requires adding one entry here.
.clusterSelectorTooltips <- c(
    selectMode = paste(
        "Choose how plot selections change the currently selected SOM nodes:",
        "view = replace; add = union; remove = subtract;",
        "remove others = intersect."
    ),
    samples2plot = paste(
        "Restrict scatter plots and statistics to the selected samples."
    ),
    scatterPercentile = paste(
        "Percentile range used to auto-zoom the scatter plot axes around the",
        "selected cells."
    ),
    somColorVar = "Statistic used to colour SOM nodes in the 2D plots.",
    somSizeVar = "Statistic used to scale SOM node points in the 2D plots.",
    clusterNumbers = paste(
        "Comma-separated SOM node ids to add to the current selection."
    ),
    applyclusterNumbers = paste(
        "Create or update the selection from the node ids above."
    ),
    clusterName = "Name for the currently selected group of nodes.",
    applyName = "Save the current selection under the name above.",
    groupRM = "Existing named groups to delete.",
    rmGroups = "Delete the selected named groups.",
    clusterNameSelect = "Existing groups to activate or view together.",
    clusterNameRM = "Single existing group to remove.",
    rmGrp = "Remove the named group selected above.",
    downloadPlots = "Download the current plots as a PDF.",
    close = "Close the clusterSelector window.",
    currentDimX = "Marker/channel for the X axis of the main SOM 2D plot.",
    currentDimY = "Marker/channel for the Y axis of the main SOM 2D plot.",
    dimPairSelect = "Preset marker pair; choosing it updates both axes.",
    showGroups = paste(
        "Colour nodes by saved group membership instead of by statistic."
    ),
    nnAgreement = paste(
        "Colour nodes by how well this marker pair preserves each node's",
        "nearest neighbours in the full SOM code space (kNN agreement)."
    ),
    colorbyGroups = "Saved groups used for node colouring.",
    dimRedSelection = "Markers used for t-SNE, UMAP and PCA calculations.",
    perplexity = "t-SNE perplexity; lower values focus on local structure.",
    n_neighbors = "UMAP number of neighbours; affects local/global balance.",
    showlegend = "Show a legend on the dimension-reduction plots.",
    staticSomX = "Marker for the X axis of a new static SOM pair view.",
    staticSomY = "Marker for the Y axis of a new static SOM pair view.",
    staticSomAdd = "Add the chosen pair to the static grid.",
    staticSomRemove = "Existing preset pair to remove from the grid.",
    staticSomRemoveBtn = "Delete the selected preset pair.",
    staticSomCols = "Number of plots per row in the static grid.",
    staticSomHeight = "Height in pixels for each static SOM plot.",
    flowSOMPieCols = "Number of pie-chart columns in the FlowSOM view.",
    flowSOMPieSize = "Pixel size of each FlowSOM pie.",
    flowSOMPieMax = "Maximum number of pies to display in the FlowSOM view.",
    compareStatsTo = paste(
        "Numerical column or saved group to compare cell counts against."
    ),
    relativeTo = "Reference population used for percentage calculations.",
    singleNode = paste(
        "Show counts for one specific SOM node alongside the selection."
    ),
    groupsVar = "Sample-level factor used to split samples for the t-test.",
    group1 = "Samples belonging to group 1 for the t-test.",
    group2 = "Samples belonging to group 2 for the t-test.",
    violinSelection = "Markers displayed in the violin plots.",
    upsetSelection = "Saved groups compared in the UpSet plot.",
    ridgeMarkerSelection = paste(
        "Markers shown as ridge rows in the marker ridge plot."
    ),
    ridgeThreshold = paste(
        "Fraction of the selected cells inside the band on each marker",
        "(0.99: the band spans the 0.5% to 99.5% quantiles); the remaining",
        "cells (1%) are shown as outlier dots."
    ),
    ridgeMaxLines = paste(
        "Maximum number of outlier cells shown; the most extreme cells",
        "(largest |z| within the selection) are kept."
    ),
    ridgeSortRows = paste(
        "Order the marker rows by how far the selection mean is shifted from",
        "the all-cell mean (in all-cell sd)."
    ),
    ridgeClusterSelection = paste(
        "When more than five clusters are selected, choose which clusters get",
        "individual density lines; the combined line always covers all",
        "selected clusters."
    ),
    diagGroups = paste(
        "Restrict the sample diagnostics plots to cells from these",
        "saved clusters/groups; leave every group selected to use the",
        "whole dataset."
    ),
    diagColorBy = "Colour the sample-level diagnostic plots by this column.",
    orderedHeatmapGroups = paste(
        "Saved clusters/groups to show as rows of the ordered heatmap."
    ),
    orderedHeatmapMarkers = paste(
        "Markers to show as columns of the ordered heatmap."
    ),
    diffGroups = paste(
        "Saved clusters/groups compared by the differential test (at",
        "least 2 required)."
    ),
    diffAnalysisType = paste(
        "Differential abundance (cell counts per group) or differential",
        "state (marker expression per group) test."
    ),
    diffFixedCol = paste(
        "Two-level experiment_info factor used as the fixed effect",
        "(the comparison of interest)."
    ),
    diffRandomCol = paste(
        "Optional experiment_info factor added as a GLMM/LMM random",
        "intercept."
    ),
    diffMarkers = "Markers tested for differential state (DS); ignored for DA.",
    runDiffTest = "Fit the diffcyt GLMM/LMM test; this can take a few seconds."
)


#' Wrap a Control in a Tooltip Div
#'
#' @param id Input id; its text in \code{.clusterSelectorTooltips} becomes
#'   the \code{title} attribute.
#' @param tag The control.
#' @param ... Further attributes of the div (e.g. \code{style}).
#'
#' @return A \code{div} tag.
#'
#' @keywords internal
.withTip <- function(id, tag, ...) {
    shiny::div(..., title = .clusterSelectorTooltips[[id]], tag)
}

#' A Column Holding One Tooltip-Wrapped Control
#'
#' @param width Column width.
#' @inheritParams .withTip
#'
#' @return A \code{column} tag.
#'
#' @keywords internal
.tipColumn <- function(width, id, tag, ...) {
    shiny::column(width = width, .withTip(id, tag, ...))
}

#' Collapsed Primary Dashboard Box
#'
#' @param title Box title.
#' @param ... Box content.
#' @param id Optional box id.
#' @param width Box width.
#' @param infoId Key into \code{.boxInfoRegistry} (see
#'   \code{R/interface_info.R}); when supplied, an info button is added next
#'   to the title that opens a modal explaining the box.
#'
#' @return A \code{shinydashboardPlus::box}.
#'
#' @keywords internal
.primaryBox <- function(title, ..., id = NULL, width = 12, infoId = NULL) {
    boxTitle <- if (is.null(infoId)) {
        title
    } else {
        shiny::tagList(title, .infoButton(infoId))
    }
    shinydashboardPlus::box(
        id = id, title = boxTitle, solidHeader = TRUE, width = width,
        status = "primary", collapsible = TRUE, collapsed = TRUE, ...
    )
}

#' Grey Italic Usage Hint
#'
#' @param ... Text pieces.
#'
#' @return A \code{p} tag.
#'
#' @keywords internal
.selectionHint <- function(...) {
    htmltools::tags$p(...,
        style = "color:#888; font-style:italic; font-size:12px; margin:4px 0;"
    )
}

#' Resizable Output
#'
#' @param out An output tag.
#'
#' @return The output wrapped by \code{shinyjqui::jqui_resizable()}.
#'
#' @keywords internal
.resizable <- function(out) shinyjqui::jqui_resizable(out)

#' Build clusterSelector Sidebar UI
#'
#' Constructs the left-hand control panel: selection mode, sample filter,
#' scatter-plot auto-zoom, SOM 2D plot aesthetics, and group management.
#'
#' Axis controls have moved to the SOM 2D plots body box to avoid the
#' selectize-dropdown overflow that occurs inside the narrow sidebar.
#'
#' @param sce Full \code{SingleCellExperiment} for sample-id choices.
#' @param outputList Named list of cluster groupings.
#' @param nPlots Kept for API compatibility; no longer controls dim-UI layout.
#'
#' @return A \code{shinydashboard::dashboardSidebar} UI object.
#'
#' @keywords internal
.buildClusterSelectorSidebar <- function(sce, outputList, nPlots) {
    shinydashboard::dashboardSidebar(
        .sidebarAssets(),
        shinydashboard::sidebarMenu(
            id = "tabs",
            shinydashboard::menuItem("Parameters", tabName = "parameters",
                icon = shiny::icon("sliders"))
        ),
        shiny::hr(),
        shiny::conditionalPanel(
            "input.tabs == 'parameters'",
            shiny::div(
                style = "padding: 0 12px;",
                .sidebarDisplayControls(sce),
                .sidebarGroupControls(outputList)
            )
        )
    )
}

#' Scripts and Styles Loaded With the Sidebar
#'
#' @return A tag list: shinyjs, the close-window extension, the pinned
#'   plot-output height and the window-resize handler.
#'
#' @keywords internal
.sidebarAssets <- function() {
    shiny::tagList(
        shinyjs::useShinyjs(),
        shinyjs::extendShinyjs(
            text = "shinyjs.closeWindow = function() { window.close(); }",
            functions = c("closeWindow")
        ),
        # shinyjqui wraps each resizable output in a .jqui-wrapper sized 100%
        # of its parent at wrap time; for outputs hidden at that moment
        # (collapsed boxes) the parent height is auto, so the wrapper becomes
        # content-driven and each render re-sizes it - an endless loop.
        # Pinning the default plot-output height in px breaks that chain.
        shiny::tags$style(shiny::HTML(
            ".tab-content .shiny-plot-output { height: 450px; }"
        )),
        # Reset widths to fill the container but never touch heights (see
        # above); heights are pinned in px by CSS instead.
        shiny::tags$script(shiny::HTML(paste(
            "(function() {",
            "  var resizeTimer = null;",
            "  function resetAndResizePlots() {",
            "    $('.ui-resizable').css({ width: '100%' });",
            "    $('.js-plotly-plot').each(function() {",
            "      Plotly.Plots.resize(this);",
            "    });",
            "  }",
            "  $(window).on('resize', function() {",
            "    clearTimeout(resizeTimer);",
            "    resizeTimer = setTimeout(resetAndResizePlots, 200);",
            "  });",
            "})();",
            sep = "\n"
        )))
    )
}

#' Sidebar Controls for Selection Mode, Samples and SOM Plot Aesthetics
#'
#' @param sce Full \code{SingleCellExperiment} for sample-id choices.
#'
#' @return A tag list.
#'
#' @keywords internal
.sidebarDisplayControls <- function(sce) {
    samples <- unique(as.character(sce$sample_id))
    stats <- c("n", "mean", "median", "rdQu", "max")
    shiny::tagList(
        .withTip("selectMode", shiny::radioButtons("selectMode",
            "Select mode", choices = c("view", "remove others", "add",
                "remove"), selected = "view", inline = FALSE)),
        .withTip("samples2plot", shiny::selectInput("samples2plot",
            "Samples to plot", choices = samples, selected = samples,
            multiple = TRUE, selectize = TRUE)),
        .withTip("scatterPercentile", shiny::sliderInput("scatterPercentile",
            "Scatter auto-zoom percentile", min = 0.5, max = 1,
            value = 0.99, step = 0.01)),
        shiny::helpText("Lower percentiles zoom in further on dense regions."),
        .withTip("somColorVar", shiny::selectInput("somColorVar",
            "SOM 2D color by", choices = stats, selected = "n",
            multiple = FALSE)),
        .withTip("somSizeVar", shiny::selectInput("somSizeVar",
            "SOM 2D size by", choices = stats, selected = "max",
            multiple = FALSE)),
        shiny::helpText(
            "Color and size summarise marker expression per SOM node."
        ),
        shiny::hr(),
        .withTip("clusterNumbers", shiny::textInput("clusterNumbers",
            "Cluster numbers", value = "1")),
        .withTip("applyclusterNumbers", shiny::actionButton(
            "applyclusterNumbers", "Apply cluster numbers")),
        shiny::hr()
    )
}

#' Sidebar Controls for Named Groups and the Window
#'
#' @param outputList Named list of cluster groupings.
#'
#' @return A tag list.
#'
#' @keywords internal
.sidebarGroupControls <- function(outputList) {
    groups <- names(outputList)
    shiny::tagList(
        .withTip("clusterName", shiny::textInput("clusterName",
            "Name selection", value = "")),
        .withTip("applyName", shiny::actionButton("applyName", "Apply name")),
        shiny::hr(),
        .withTip("groupRM", shiny::selectInput("groupRM", "Select to remove",
            choices = groups, multiple = TRUE, selectize = TRUE)),
        .withTip("rmGroups", shiny::actionButton("rmGroups",
            "Remove groups")),
        shiny::hr(),
        .withTip("clusterNameSelect", shiny::selectInput("clusterNameSelect",
            "Select named", choices = groups, multiple = TRUE,
            selectize = TRUE)),
        .withTip("clusterNameRM", shiny::selectInput("clusterNameRM",
            "Remove by name", choices = groups, multiple = FALSE,
            selectize = TRUE)),
        .withTip("rmGrp", shiny::actionButton("rmGrp", "Remove")),
        shiny::hr(),
        .withTip("downloadPlots", shiny::downloadButton("downloadPlots",
            "Download Plots")),
        .withTip("close", shiny::actionButton("close", "Close window"))
    )
}


#' Build the First Dashboard Body Row
#'
#' Top row of the dashboard: collapsible sample tree (optional), the main 2D
#' scatter plot, and the interactive dendrogram. All three outputs support
#' node/point selection.
#'
#' @param colTree Optional collapsible tree object.
#' @param nPlots Number of 2D SOM plots (kept for signature consistency).
#'
#' @return A \code{shiny::fluidRow} object.
#'
#' @keywords internal
.buildFirstBodyRow <- function(colTree, nPlots) {
    hint <- .selectionHint(
        "Tip: click or lasso points/nodes to select clusters.",
        "Selection mode is",
        "controlled from the sidebar."
    )
    shiny::fluidRow(
        if (!is.null(colTree)) {
            shiny::column(width = 3, .primaryBox("interactive tree",
                collapsibleTree::collapsibleTreeOutput("plot"),
                infoId = "interactiveTree"))
        },
        shiny::column(width = 6, .primaryBox("2D plot", hint,
            .resizable(plotly::plotlyOutput("scatter")),
            infoId = "twoDPlot")),
        shiny::column(width = 6, .primaryBox("interactive dendrogram", hint,
            .resizable(plotly::plotlyOutput("dendPlotly")),
            infoId = "dendrogram"))
    )
}

#' Initial Axis Pair
#'
#' @param markers Available markers.
#' @param dList List of marker pairs.
#'
#' @return The first preset pair, or the first two markers.
#'
#' @keywords internal
.initialAxes <- function(markers, dList) {
    if (length(dList) > 0L) {
        return(dList[[1L]][seq_len(2L)])
    }
    c(markers[1L], markers[min(2L, length(markers))])
}

#' Build SOM 2D Plot Box
#'
#' Main SOM visualisation panel: interactive axis controls, the primary 2D SOM
#' plot, dimension-reduction controls, and t-SNE/UMAP/PCA projection plots.
#' All plot outputs are selection-enabled.
#'
#' @param nPlots Kept for API compatibility.
#' @param colsUsed Character vector of SOM column names for dim-red selector.
#' @param outputList Named list of cluster groupings.
#' @param markers Character vector of all available markers
#'   (\code{rownames(sce)}).
#' @param dList List of marker pairs used to populate the preset-pair picker.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildSOM2DPlotsBox <- function(nPlots, colsUsed, outputList, markers, dList) {
    init <- .initialAxes(markers, dList)
    hint <- .selectionHint(
        "Tip: click or lasso points to select SOM nodes. Selection mode is",
        "controlled from the sidebar."
    )
    axisSelect <- function(id, label, selected) {
        .tipColumn(4L, id, shiny::selectInput(id, label, choices = markers,
            selected = selected, multiple = FALSE, selectize = TRUE,
            width = "100%"))
    }
    .primaryBox(
        "SOM 2D plots", width = 12L, infoId = "som2dPlots",
        shiny::fluidRow(
            axisSelect("currentDimX", "X axis", init[1L]),
            axisSelect("currentDimY", "Y axis", init[2L]),
            ## Rendered by server so it stays in sync with dListRV
            .tipColumn(4L, "dimPairSelect", shiny::uiOutput("dimPairSelectUI"))
        ),
        shiny::fluidRow(
            .tipColumn(2L, "showGroups",
                shiny::checkboxInput("showGroups", "Colour groups",
                    value = FALSE)),
            .tipColumn(2L, "nnAgreement",
                shiny::checkboxInput("nnAgreement",
                    "Neighbourhood agreement", value = FALSE))
        ),
        shiny::fluidRow(shiny::column(width = 12L, hint,
            .resizable(plotly::plotlyOutput("somDataMain")))),
        .dimRedControls(colsUsed, outputList),
        shiny::fluidRow(lapply(c("tsne", "umap", "pca"), function(id) {
            shiny::column(width = 4L, hint,
                .resizable(plotly::plotlyOutput(id)))
        }))
    )
}

#' Dimension-Reduction Controls of the SOM 2D Plot Box
#'
#' @param colsUsed Character vector of SOM column names.
#' @param outputList Named list of cluster groupings.
#'
#' @return A \code{fluidRow}.
#'
#' @keywords internal
.dimRedControls <- function(colsUsed, outputList) {
    cn <- setdiff(colsUsed, c("label", "clusterid"))
    shiny::fluidRow(
        shiny::column(
            width = 8L,
            .withTip("colorbyGroups", shiny::selectInput("colorbyGroups",
                "Colour by groups", choices = names(outputList),
                multiple = TRUE, selectize = TRUE)),
            .withTip("dimRedSelection", shiny::selectInput("dimRedSelection",
                "Markers for dim. reduction", choices = cn, selected = cn,
                multiple = TRUE, selectize = TRUE)),
            shiny::helpText(paste("Dimension-reduction runs on SOM code",
                "vectors using the selected markers."))
        ),
        shiny::column(
            width = 2L,
            .withTip("perplexity", shiny::numericInput("perplexity",
                "Perplexity", value = 30, min = 1, max = 500)),
            .withTip("showlegend", shiny::checkboxInput("showlegend",
                "Show legend", value = FALSE))
        ),
        .tipColumn(2L, "n_neighbors", shiny::numericInput("n_neighbors",
            "n_neighbors", value = 4, min = 2, max = 500))
    )
}


#' Build Six Static SOM Views Box
#'
#' Grid of static ggplot SOM pair views. Users can add/remove marker pairs,
#' change the grid layout, and adjust per-plot height.
#'
#' @param markers Character vector of all available markers
#'   (\code{rownames(sce)}).
#' @param dList Initial list of marker pairs.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#' @keywords internal
.buildSixStaticSOMBox <- function(markers, dList) {
    .primaryBox(
        "All SOM pair views (static)", width = 12L,
        infoId = "staticSomViews",
        .staticSomPairControls(markers, dList),
        shiny::fluidRow(
            .tipColumn(3L, "staticSomCols", shiny::selectInput(
                "staticSomCols", "Columns per row",
                choices = c("1" = 1L, "2" = 2L, "3" = 3L, "4" = 4L, "6" = 6L),
                selected = 3L, multiple = FALSE, selectize = FALSE,
                width = "100%")),
            .tipColumn(3L, "staticSomHeight", shiny::sliderInput(
                "staticSomHeight", "Plot height (px)",
                min = 150L, max = 600L, value = 240L, step = 10L))
        ),
        shiny::hr(),
        htmltools::tags$p(
            "Red dots = selected SOM nodes. Expand box to render.",
            style = paste(
                "color:#888; font-style:italic;",
                "font-size:12px; margin-bottom:8px;"
            )
        ),
        ## Dynamic plot grid (filled by renderUI in server).
        shiny::uiOutput("staticSomGrid")
    )
}

#' Add / Remove Controls of the Static SOM Pair Grid
#'
#' @param markers Available markers.
#' @param dList Initial list of marker pairs.
#'
#' @return A \code{fluidRow}.
#'
#' @keywords internal
.staticSomPairControls <- function(markers, dList) {
    removeChoices <- if (length(dList) > 0L) {
        stats::setNames(seq_along(dList), vapply(dList,
            function(p) paste(p[1L], "-", p[2L]), character(1L)))
    } else {
        character(0L)
    }
    init <- .initialAxes(markers, dList)
    axisSelect <- function(id, label, selected) {
        .tipColumn(3L, id, shiny::selectInput(id, label, choices = markers,
            selected = selected, multiple = FALSE, selectize = TRUE,
            width = "100%"))
    }
    button <- function(width, id, label, icon, class) {
        .tipColumn(width, id, shiny::actionButton(id, label = label,
            icon = shiny::icon(icon), class = class),
            style = "margin-top:25px;")
    }
    shiny::fluidRow(
        axisSelect("staticSomX", "X axis (add)", init[1L]),
        axisSelect("staticSomY", "Y axis (add)", init[2L]),
        button(2L, "staticSomAdd", "Add", "plus", "btn-success btn-sm"),
        # Native <select> (selectize = FALSE) never overflows the box.
        .tipColumn(3L, "staticSomRemove", shiny::selectInput(
            "staticSomRemove", "Remove pair", choices = removeChoices,
            multiple = FALSE, selectize = FALSE, width = "100%")),
        button(1L, "staticSomRemoveBtn", NULL, "trash", "btn-danger btn-sm")
    )
}

#' Build FlowSOM Marker Pies UI
#'
#' Displays each selected SOM node as a pie chart coloured by mean marker
#' expression. Users control the grid layout and pie size.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#' @keywords internal
.buildFlowSOMPieBox <- function() {
    .primaryBox(
        "FlowSOM marker pies", width = 12L, infoId = "flowSomPies",
        shiny::fluidRow(
            .tipColumn(3L, "flowSOMPieCols", shiny::numericInput(
                inputId = "flowSOMPieCols", label = "Columns",
                value = 5L, min = 1L, max = 20L, step = 1L, width = "100%")),
            .tipColumn(3L, "flowSOMPieSize", shiny::sliderInput(
                inputId = "flowSOMPieSize", label = "Size per pie (px)",
                min = 50L, max = 500L, value = 150L, step = 10L,
                width = "100%")),
            .tipColumn(3L, "flowSOMPieMax", shiny::numericInput(
                inputId = "flowSOMPieMax", label = "Max pies",
                value = 5L, min = 1L, max = 50L, step = 1L, width = "100%")),
            shiny::column(width = 3L, htmltools::tags$p(
                paste(
                    "Each pie = one selected SOM node.",
                    "Slice size = mean marker expression.",
                    "Colours are consistent across nodes."
                ),
                style = "color:#666; font-size:12px; margin-top:8px;"
            ))
        ),
        ## Height is set dynamically in the server so each pie keeps a
        ## consistent size regardless of how many nodes are selected.
        shiny::fluidRow(shiny::column(width = 12L,
            shiny::uiOutput("flowSOMPieUI")))
    )
}

#' Build FlowSOM Stars UI
#'
#' Interactive FlowSOM star plot. Each node shows marker expression as a star,
#' background colour shows the metaclustering. Click or lasso nodes to select
#' clusters.
#'
#' @param fsom Optional \code{FlowSOM} object. If \code{NULL}, returns
#'   \code{NULL}.
#'
#' @return A \code{shinydashboardPlus::box} UI object or \code{NULL}.
#' @keywords internal
.buildFlowSOMStarsBox <- function(fsom) {
    if (is.null(fsom)) {
        return(NULL)
    }
    .primaryBox(
        "FlowSOM stars", infoId = "flowSomStars",
        shiny::fluidRow(shiny::column(
            width = 12,
            htmltools::tags$p(
                "Interactive FlowSOM star plot. Each node shows marker",
                "expression",
                "as a star; background colour shows metaclustering. Click or",
                "lasso nodes to select clusters."
            ),
            .resizable(plotly::plotlyOutput("flowSOMStars"))
        ))
    )
}


#' Build Stats UI
#'
#' Statistics panel: compare cell counts, compute percentages, run sample-group
#' t-tests, and display the results as tables and bar plots.
#'
#' @param metaD \code{metadata(sce)} list used to build numeric/factor choices.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param outputList Named list of cluster groupings.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildStatsBox <- function(metaD, somCodesName, outputList) {
    info <- metaD$experiment_info
    numCols <- unlist(lapply(info, is.numeric), use.names = FALSE)
    eI <- apply(info[, numCols, drop = FALSE], 2, as.numeric)
    if (any(is.na(eI))) {
        stop("NAs produced when converting experiment_info to numeric")
    }
    factCols <- unlist(lapply(info, is.factor), use.names = FALSE)
    .primaryBox(
        "Stats", infoId = "stats",
        shiny::fluidRow(
            shiny::column(width = 4, shiny::verbatimTextOutput("somClusters")),
            shiny::column(
                width = 8,
                .statsControls(c("none", colnames(eI)),
                    c("none", colnames(info)[factCols]),
                    nrow(metaD[[somCodesName]])),
                .statsResultBoxes()
            )
        )
    )
}

#' Comparison and t-Test Controls of the Stats Box
#'
#' @param numChoices Numeric \code{experiment_info} columns (plus "none").
#' @param factChoices Factor \code{experiment_info} columns (plus "none").
#' @param nNodes Number of SOM nodes.
#'
#' @return A tag list of two fluid rows.
#'
#' @keywords internal
.statsControls <- function(numChoices, factChoices, nNodes) {
    pick <- function(id, label, choices, selected, multiple) {
        .tipColumn(4, id, shiny::selectInput(id, label, choices = choices,
            selected = selected, multiple = multiple, selectize = TRUE))
    }
    shiny::tagList(
        shiny::fluidRow(
            pick("compareStatsTo", "select comparison Stats", numChoices,
                "none", FALSE),
            pick("relativeTo", "Parent population for percentages",
                numChoices, "none", FALSE),
            .tipColumn(4, "singleNode", shiny::numericInput("singleNode",
                "single SOM node counts", min = 1, max = nNodes, value = 1,
                step = 1))
        ),
        shiny::fluidRow(
            pick("groupsVar", "Select groups variable for t-test",
                factChoices, "none", FALSE),
            pick("group1", "group 1", c(), "", TRUE),
            pick("group2", "group 2", c(), "", TRUE)
        )
    )
}

#' Count Table, Bar Plots and t-Test Result Boxes of the Stats Box
#'
#' @return A tag list of two fluid rows.
#'
#' @keywords internal
.statsResultBoxes <- function() {
    shiny::tagList(
        shiny::fluidRow(.primaryBox(
            "Stats",
            shiny::fluidRow(
                .resizable(DT::DTOutput("cellCounts")),
                .resizable(shiny::verbatimTextOutput("cellPercentages"))
            ),
            shiny::fluidRow(
                shiny::column(width = 6,
                    .resizable(shiny::plotOutput("CountBar"))),
                shiny::column(width = 6,
                    .resizable(shiny::plotOutput("PercentBar")))
            )
        )),
        shiny::fluidRow(.primaryBox(
            "t-Test", infoId = "tTest",
            shiny::fluidRow(shiny::column(width = 12,
                shiny::verbatimTextOutput("ttestResult")))
        ))
    )
}


#' Build SOM Raster UI
#'
#' Raster heatmap of SOM node marker expression plus a small interactive plot
#' for selecting nodes from the SOM grid.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildSOMRasterBox <- function() {
    .primaryBox(
        "som plots", infoId = "somRaster",
        shiny::fluidRow(shiny::column(width = 12,
            .resizable(shiny::plotOutput("somRaster", height = "1200px")))),
        shiny::fluidRow(shiny::column(
            width = 12,
            .selectionHint(
                "Tip: click or lasso nodes in the right-hand plot to add them",
                "to the current selection."
            ),
            .resizable(plotly::plotlyOutput("somRasterSelect",
                height = "400px"))
        ))
    )
}


#' Build Violin Plots UI
#'
#' Marker expression violin plots for the saved cluster groups. Users choose
#' which markers are shown.
#'
#' @param somCodes SOM codes matrix; column names become marker choices.
#' @param colsUsed Default marker selection.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildViolinBox <- function(somCodes, colsUsed) {
    .primaryBox(
        "violin plots", id = "violinBox", infoId = "violinPlots",
        shiny::fluidRow(
            shiny::column(width = 12),
            .withTip("violinSelection", shiny::selectInput("violinSelection",
                "Select markers to show", choices = colnames(somCodes),
                selected = colsUsed, multiple = TRUE, selectize = TRUE))
        ),
        shiny::fluidRow(
            shiny::column(width = 6, .resizable(shiny::plotOutput("VlnPlot"))),
            shiny::column(width = 6, .resizable(shiny::plotOutput("VlnPlot2")))
        )
    )
}


#' Build UpSet Plot UI
#'
#' UpSet-style overlap visualisation of saved cluster groups. Users choose which
#' groups to compare.
#'
#' @param outputList Named list of cluster groupings.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildUpSetBox <- function(outputList) {
    default <- intersect(c("Rest", "selected"), names(outputList))
    if (length(default) == 0) default <- names(outputList)
    .primaryBox(
        "UpSet plot", id = "upsetBox", infoId = "upsetPlot",
        shiny::fluidRow(
            shiny::column(width = 12),
            .withTip("upsetSelection", shiny::selectInput("upsetSelection",
                "Select groups to show", choices = names(outputList),
                selected = default, multiple = TRUE, selectize = TRUE))
        ),
        shiny::fluidRow(shiny::column(width = 12,
            .resizable(shiny::plotOutput("UpSet"))))
    )
}


#' Build Marker Ridge Trace Box UI
#'
#' Interactive per-marker ridgeline densities of all cells and of the
#' current selection, with a band holding the central fraction (coverage)
#' of the selected cells and the selected cells outside that band as
#' clickable outlier dots; a second plot explains the SOM node membership of
#' the clicked cell. Users choose the markers, the band coverage and the
#' outlier-cell cap.
#'
#' @param somCodes SOM codes matrix; column names become marker choices.
#' @param colsUsed Default marker selection.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildRidgeBox <- function(somCodes, colsUsed) {
    .primaryBox(
        "marker ridge", id = "ridgeBox", infoId = "ridgeBox",
        shiny::fluidRow(
            .tipColumn(4, "ridgeMarkerSelection", shiny::selectInput(
                "ridgeMarkerSelection", "Select markers to show",
                choices = colnames(somCodes), selected = colsUsed,
                multiple = TRUE, selectize = TRUE)),
            .tipColumn(4, "ridgeThreshold", shiny::numericInput(
                "ridgeThreshold", "Band coverage (fraction of selected cells)",
                value = 0.99, min = 0.5, max = 1, step = 0.005)),
            .tipColumn(4, "ridgeMaxLines", shiny::numericInput(
                "ridgeMaxLines", "Max outlier cells",
                value = 50, min = 1, max = 500, step = 1))
        ),
        shiny::fluidRow(
            # Cluster picker: filled by the server when more than five
            # clusters are selected (individual lines are capped at five).
            shiny::column(width = 8, shiny::uiOutput("ridgeClusterUI")),
            .tipColumn(4, "ridgeSortRows", shiny::checkboxInput(
                "ridgeSortRows", "Sort rows by separation", value = TRUE))
        ),
        shiny::fluidRow(shiny::column(width = 12, .resizable(
            plotly::plotlyOutput("ridgePlot", height = "600px")))),
        # Why does the clicked cell belong to its SOM node?
        shiny::fluidRow(shiny::column(width = 12, .resizable(
            plotly::plotlyOutput("ridgeCellWhy", height = "300px"))))
    )
}
