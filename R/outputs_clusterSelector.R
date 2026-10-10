# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# outputs_clusterSelector.R ----
# Shiny output renderers for the clusterSelector app. Outputs are registered
# separately from observers to follow Bioconductor's Shiny app packaging
# guidelines. Renderers read the per-session state environment `st` built in
# server_clusterSelector.R; the plot and table logic lives in the helpers of
# outputs-helpers.R.

#' Register clusterSelector Outputs
#'
#' Attaches all `output$*` renderers used by the clusterSelector app. Observer
#' logic lives in `.registerClusterSelectorObservers()`.
#'
#' @param input,output,session Shiny objects.
#' @param st Server state environment (see `.initServerState()`).
#'
#' @return Called for its side effects.
#'
#' @keywords internal
#' @name INTERNAL_registerClusterSelectorOutputs
.registerClusterSelectorOutputs <- function(input, output, session, st) {
    .registerSelectionOutputs(input, output, st)
    .registerSOMOutputs(input, output, st)
    .registerStaticSOMOutputs(input, output, st)
    .registerProjectionOutputs(input, output, st)
    .registerStatsOutputs(input, output, st)
    .registerRasterPieOutputs(input, output, st)
    .registerGroupPlotOutputs(input, output, st)
    .registerDiagnosticsOutputs(input, output, st)
    .registerOrderedHeatmapOutputs(input, output, st)
    .registerDiffTestOutputs(input, output, st)
    output$downloadPlots <- .downloadPlotsHandler(input, st)
    invisible(NULL)
}

#' Register the Sample Tree, Preset Pair Picker and Cell Scatter Outputs
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerSelectionOutputs <- function(input, output, st) {
    if (!is.null(st$colTree)) {
        output$plot <- collapsibleTree::renderCollapsibleTree(st$colTree)
    }
    # Dimension-pair preset picker rendered from the reactive dListRV.
    output$dimPairSelectUI <- shiny::renderUI({
        dList <- st$dListRV()
        if (length(dList) == 0L) {
            return(shiny::selectInput("dimPairSelect", "Preset pair",
                choices = character(0L)))
        }
        labels <- vapply(dList, function(p) paste(p[1L], "-", p[2L]),
            character(1L))
        shiny::selectInput(
            "dimPairSelect", "Preset pair",
            choices = stats::setNames(seq_along(dList), labels),
            selected = shiny::isolate(st$activePlot())
        )
    })
    output$scatter <- plotly::renderPlotly(.scatterPlotly(st, input))
    invisible(NULL)
}

#' Register the SOM 2D Plot Outputs
#'
#' One plotly output per SOM plot (\code{somData1}, ...) plus the single
#' main plot (\code{somDataMain}) of the active axis pair.
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerSOMOutputs <- function(input, output, st) {
    lapply(seq_len(st$nPlots), function(i) {
        output[[paste0("somData", i)]] <- plotly::renderPlotly(
            .somPlotly(st, input, i)
        )
    })
    output$somDataMain <- plotly::renderPlotly(.somMainPlot(st, input))
    invisible(NULL)
}

#' Register the Static SOM Pair Grid Outputs
#'
#' Renders the grid layout from \code{dListRV} and one static plot per
#' preset pair (\code{staticSomDyn<i>}).
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerStaticSOMOutputs <- function(input, output, st) {
    output$staticSomGrid <- shiny::renderUI(
        .staticSomGridUI(st$dListRV(), input$staticSomCols,
            input$staticSomHeight)
    )
    shiny::observe({
        current <- st$dListRV()
        lapply(seq_along(current), function(i) {
            outId <- paste0("staticSomDyn", i)
            output[[outId]] <- shiny::renderPlot(
                .staticSomPlot(st, input, current[[i]])
            )
            shiny::outputOptions(output, outId, suspendWhenHidden = TRUE)
        })
    })
    invisible(NULL)
}

#' Register the Dimension-Reduction and Dendrogram Outputs
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerProjectionOutputs <- function(input, output, st) {
    output$tsne <- plotly::renderPlotly(
        .projectionPlotly(st$tsnePlot(), "tsne", input$showlegend)
    )
    output$umap <- plotly::renderPlotly(
        .projectionPlotly(st$umapPlot(), "umap", input$showlegend)
    )
    output$pca <- plotly::renderPlotly(
        .projectionPlotly(st$pcaPlot(), "pca", input$showlegend)
    )
    output$dend <- shiny::renderPlot({
        st$dendPlot() %>% plot(main = "dendrogram")
    })
    output$dendPlotly <- plotly::renderPlotly({
        dd <- st$dendPlotlyData()
        rs <- st$rsUsed()
        shiny::req(dd, rs)
        .dendPlotly(dd, rs)
    })
    invisible(NULL)
}

#' Register the Statistics Outputs
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerStatsOutputs <- function(input, output, st) {
    output$somClusters <- shiny::renderText({
        rs <- st$rsUsed()
        shiny::req(rs)
        .formatNodeList(as.integer(
            intersect(rs, colnames(st$clusterPatientTable))
        ))
    })
    output$cellCounts <- DT::renderDT(
        options = list(lengthChange = FALSE, scrollX = TRUE),
        .cellCountsTable(st, input)
    )
    output$cellPercentages <- shiny::renderPrint(
        .cellPercentages(st, input)
    )
    output$CountBar <- shiny::renderPlot(st$countBarPlot())
    output$PercentBar <- shiny::renderPlot(st$PercentBarPlot())
    output$ttestResult <- shiny::renderPrint(.ttestResult(st, input))
    invisible(NULL)
}

#' Register the SOM Raster and FlowSOM Outputs
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerRasterPieOutputs <- function(input, output, st) {
    output$somRaster <- shiny::renderPlot({
        xy <- st$somRasterPlot()
        shiny::req(xy, st$baseRasterGgplot)
        .rasterWithSelection(st$baseRasterGgplot, xy)
    })
    output$somRasterSelect <- plotly::renderPlotly({
        rs <- st$rsUsed()
        if (is.null(rs)) {
            return(NULL)
        }
        p3 <- .buildSOMRasterSelectPlot(st$somRasterData, rs)
        shiny::req(inherits(p3, "ggplot"))
        quiet_ggplotly(p3, source = "somGrid", tooltip = "") %>%
            plotly::layout(showlegend = FALSE, dragmode = "select") %>%
            plotly::event_register("plotly_selected") %>%
            plotly::event_register("plotly_relayout")
    })
    output$flowSOMPieUI <- shiny::renderUI({
        rs <- st$rsUsed()
        shiny::req(rs)
        # Size the output for the pies actually drawn (at most Max pies).
        nPies <- min(length(rs), input$flowSOMPieMax %||% 5L)
        .flowSOMPieUI(nPies, input$flowSOMPieCols, input$flowSOMPieSize)
    })
    if (!is.null(st$fsom)) {
        output$flowSOMStars <- plotly::renderPlotly(
            .flowSOMStarsPlotly(st$fsom)
        )
    }
    output$flowSOMPie <- shiny::renderPlot({
        p <- st$flowSOMPiePlot()
        shiny::req(p)
        p
    })
    invisible(NULL)
}

#' Register the Violin, UpSet and Ridge Cluster-Picker Outputs
#'
#' @param input,output Shiny objects.
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerGroupPlotOutputs <- function(input, output, st) {
    output$VlnPlot <- shiny::renderPlot(st$vlnPlot())
    output$VlnPlot2 <- shiny::renderPlot(st$VlnPlot2())
    output$UpSet <- shiny::renderPlot(st$upSetPlot())
    # Cluster picker for the marker ridge box: shown only when more than
    # five clusters are selected, because individual density lines are
    # capped at five. Selecting none keeps the combined line only.
    output$ridgeClusterUI <- shiny::renderUI({
        rs <- st$rsUsed()
        shiny::req(rs)
        .ridgeClusterUI(st$sce, rs, input$ridgeClusterSelection)
    })
    invisible(NULL)
}

#' PDF Download Handler for All Plots
#'
#' @param input Shiny input object.
#' @param st Server state environment.
#'
#' @return A \code{shiny::downloadHandler}.
#'
#' @keywords internal
.downloadPlotsHandler <- function(input, st) {
    shiny::downloadHandler(
        filename = function() {
            nm <- input$clusterNameSelect
            if (is.null(nm) || !nzchar(nm)) nm <- "clusterSelector_plots"
            paste0(nm, ".pdf")
        },
        content = function(file) {
            rs <- st$rsUsed()
            shiny::req(rs)
            .writeClusterSelectorPdf(
                file = file, rs = rs, sce = st$sce, sceRN = st$sceRN,
                sceCN = st$sceCN, metaD = st$metaD,
                somCodesName = st$somCodesName,
                dendPlotObj = st$dendPlot(),
                dimPairs = st$dListRV(), getBasePlotFn = st$getBasePlot,
                dlColorVar = input$somColorVar %||% "n",
                dlSizeVar = input$somSizeVar %||% "max",
                pctl = input$scatterPercentile %||% 0.99,
                somRasterXy = st$somRasterPlot(),
                baseRasterGgplot = st$baseRasterGgplot,
                plotRegistry = .pdfPlotRegistry(st, input, rs)
            )
        }
    )
}
