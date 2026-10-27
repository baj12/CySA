# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server-reactives.R ----
# Reactive builders of the clusterSelector server. Each builder adds its
# reactives to the server state environment `st` (see
# server_clusterSelector.R); reactive bodies read other state lazily via
# `st$<name>`.

#' Build the Axis and Dimension-Reduction Reactives
#'
#' Adds \code{activeDims} (debounced axis pair), the observer keeping
#' \code{dimSelection} in sync with it, the dimension-reduction reactives
#' \code{dimRedSelection}, \code{tsne}, \code{umap}, \code{pca} and their
#' plots \code{tsnePlot}, \code{umapPlot}, \code{pcaPlot} to \code{st}.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Called for its side effect on \code{st}.
#'
#' @keywords internal
.buildDimRedPlotReactives <- function(st, input) {
    st$activeDims <- shiny::reactive({
        shiny::req(input$currentDimX, input$currentDimY)
        c(input$currentDimX, input$currentDimY)
    }) %>% shiny::debounce(250)
    shiny::observe(.syncDimSelection(st, st$activeDims()))
    codes <- function() st$metaD[[st$somCodesName]]
    st$dimRedSelection <- shiny::reactive(input$dimRedSelection) %>%
        shiny::debounce(1000)
    st$tsne <- shiny::reactive({
        tsneFunc(dimRedSelection = st$dimRedSelection(),
            perplexity = input$perplexity, st$sce, st$somCodesName)
    }) %>%
        shiny::bindCache(st$dimRedSelection(), input$perplexity) %>%
        shiny::debounce(1000)
    st$umap <- shiny::reactive({
        pumap <- umap::umap.defaults
        pumap$n_neighbors <- input$n_neighbors
        umap::umap(codes()[, input$dimRedSelection], config = pumap)
    }) %>%
        shiny::bindCache(input$dimRedSelection, input$n_neighbors) %>%
        shiny::debounce(1000)
    st$pca <- shiny::reactive({
        stats::prcomp(t(codes()[, input$dimRedSelection]),
            scale = FALSE, rank. = 2)
    }) %>%
        shiny::bindCache(input$dimRedSelection) %>%
        shiny::debounce(1000)
    st$tsnePlot <- .projectionPlotReactive(st, input,
        function() st$tsne()$Y, c("tsne1", "tsne2"))
    st$umapPlot <- .projectionPlotReactive(st, input,
        function() st$umap()$layout, c("umap1", "umap2"))
    st$pcaPlot <- .projectionPlotReactive(st, input,
        function() st$pca()$rotation, c("pc1", "pc2"))
    invisible(st)
}

#' Reset the Dimension Selection When the Axis Pair Changes
#'
#' @param st Server state environment.
#' @param dims Character vector of the two active channels.
#'
#' @return Called for its side effect on \code{st$dimSelection}.
#'
#' @keywords internal
.syncDimSelection <- function(st, dims) {
    shiny::req(dims[1L], dims[2L])
    lim1 <- st$channelLimits[[dims[1L]]]
    lim2 <- st$channelLimits[[dims[2L]]]
    shiny::req(lim1, lim2)
    existing <- shiny::isolate(st$dimSelection())
    if (length(existing) >= 1L && identical(existing[[1L]]$dims, dims)) {
        return(invisible(NULL))
    }
    st$dimSelection(list(list(
        dims = dims,
        xlim = c(lim1["min"], lim1["max"]),
        ylim = c(lim2["min"], lim2["max"]),
        xzoom = c(NULL, NULL),
        yzoom = c(NULL, NULL)
    )))
}

#' Reactive Projection Plot (t-SNE / UMAP / PCA)
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#' @param coords Function returning the two-column coordinate matrix of
#'   the SOM nodes.
#' @param axisNames Column names for the coordinates.
#'
#' @return A reactive returning the projection ggplot.
#'
#' @keywords internal
.projectionPlotReactive <- function(st, input, coords, axisNames) {
    shiny::reactive({
        xy <- coords()
        st$selectedUpdate2()
        rs <- st$rsUsed_d()
        shiny::req(rs, xy)
        st$triggerRedraw()
        df <- as.data.frame(xy)
        names(df) <- axisNames
        drawProjection(df, rs,
            colorbyGroups = input$colorbyGroups,
            sce = st$sce, outputList = st$rv$outputList
        )
    })
}

#' Build the Selection Plot Reactives
#'
#' Adds the SOM base plot and projection data, the cell scatter plot, the
#' SOM raster selection and the dendrogram reactives to \code{st}.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Called for its side effect on \code{st}.
#'
#' @keywords internal
.buildSelectionPlotReactives <- function(st, input) {
    st$somBasePlot <- shiny::reactive({
        dims <- st$activeDims()
        shiny::req(dims)
        st$getBasePlot(dims[1L], dims[2L], input$somColorVar %||% "n",
            input$somSizeVar %||% "max")
    }) %>%
        shiny::bindCache(
            st$activeDims(), input$somColorVar, input$somSizeVar
        )
    st$somProjectionDf <- shiny::reactive({
        pp1 <- st$somBasePlot()
        dimSel <- st$dimSelection()
        shiny::req(pp1, length(dimSel) >= 1L)
        buildProjectionDf(pp1, 1L, dimSel, st$sce,
            somCodesName = st$somCodesName)
    }) %>% shiny::bindCache(st$activeDims())
    # Per-node kNN agreement of the displayed pair versus the full SOM
    # code space, for the optional neighbourhood-agreement view.
    st$nnAgreement <- shiny::reactive({
        df <- st$somProjectionDf()
        shiny::req(df)
        codes <- st$metaD[[st$somCodesName]]
        xy <- as.matrix(df[, seq_len(2), drop = FALSE])
        .knnAgreement(codes, xy, k = min(10L, max(1L, nrow(codes) - 1L)))
    }) %>% shiny::bindCache(st$activeDims())
    st$scatterPlot <- shiny::reactive(.scatterPlotFromState(st, input))
    st$somRasterPlot <- shiny::reactive({
        rs <- as.integer(unlist(st$rsUsed()))
        rs <- rs[!is.na(rs) & rs > 0L]
        shiny::req(length(rs) > 0L)
        st$somRasterData[st$somRasterData$id %in% rs, c("x", "y"),
            drop = FALSE]
    })
    st$dendPlot <- shiny::reactive({
        rs <- st$rsUsed()
        shiny::req(rs)
        st$selectedPoints()
        labCol <- rep("blue", length(labels(st$dend)))
        labCol[which(labels(st$dend) %in% rs)] <- "red"
        st$dend %>%
            dendextend::set("leaves_col", labCol) %>%
            dendextend::set("leaves_pch", 15) %>%
            dendextend::set("leaves_cex", 1) %>%
            dendextend::set("labels_cex", 0.5)
    })
    st$dendPlotlyData <- shiny::reactive(
        .dendPlotlyData(st$dend, st$dendTable)
    )
    invisible(st)
}

#' Cell Scatter Plot of the Selected Nodes
#'
#' Hexbin scatter of the subsampled cells of the selected SOM nodes and
#' samples, zoomed to the \code{input$scatterPercentile} range.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return A ggplot, or \code{NULL} when no cell is selected.
#'
#' @keywords internal
.scatterPlotFromState <- function(st, input) {
    rs <- st$rsUsed()
    shiny::req(rs)
    dimSel <- st$dimSelection()
    sampleIds <- input$samples2plot
    shiny::req(sampleIds)
    st$activePlot()
    shiny::req(dimSel, length(dimSel) >= 1L)
    cd <- SingleCellExperiment::colData(st$sce_subsampled)
    cidIdx <- cd$cluster_id %in% rs & cd$sample_id %in% sampleIds
    if (sum(cidIdx) < 1) {
        return(NULL)
    }
    chs <- dimSel[[1L]]$dims
    lims <- .computeScatterLimits(
        st$df[cidIdx, make.names(chs[1]), drop = TRUE],
        st$df[cidIdx, make.names(chs[2]), drop = TRUE],
        input$scatterPercentile %||% 0.99
    )
    plotCytoScatter(
        rownms = st$sce_subsampledRN,
        x = st$sce_subsampled[, cidIdx],
        chs = chs,
        bins = 200
    ) +
        ggplot2::xlim(lims$xlim) +
        ggplot2::ylim(lims$ylim)
}

#' Dendrogram Segments and Labels for the Plotly Dendrogram
#'
#' @param dend Dendrogram object.
#' @param dendTable Optional data frame with \code{id} and \code{label}.
#'
#' @return A list with \code{segments} and \code{labels} (with node
#'   \code{id}).
#'
#' @keywords internal
.dendPlotlyData <- function(dend, dendTable) {
    g <- dendextend::as.ggdend(dend)
    labels <- g$labels
    if (!is.null(dendTable) &&
        all(c("id", "label") %in% colnames(dendTable))) {
        labels <- dplyr::left_join(
            labels, dendTable[, c("id", "label")],
            by = "label"
        )
    } else {
        labels$id <- .asIntegerQuiet(labels$label)
    }
    list(segments = g$segments, labels = labels)
}

#' Build the Statistics Reactives
#'
#' Adds the count and percentage bar plots, the debounced sample-group
#' selection (\code{groupsInput}) and the typed node ids
#' (\code{inputClusterNumber}) to \code{st}.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Called for its side effect on \code{st}.
#'
#' @keywords internal
.buildStatsReactives <- function(st, input) {
    st$countBarPlot <- shiny::reactive({
        rs <- st$rsUsed()
        shiny::req(rs, st$clusterPatientTable)
        countBarPlotFunc(rs, st$clusterPatientTable, input$compareStatsTo,
            st$sce, st$rv$outputList, st$groupsInput())
    })
    st$PercentBarPlot <- shiny::reactive({
        rs <- st$rsUsed()
        shiny::req(rs, st$clusterPatientTable)
        PercentBarPlotFunc(
            st$sce, input$relativeTo, st$clusterPatientTable, rs,
            st$rv$outputList, group, st$groupsInput()
        )
    })
    st$groupsInput <- shiny::reactive(
        .groupsFromInput(st$metaD, input$groupsVar, input$group1,
            input$group2)
    ) %>% shiny::debounce(1000)
    st$inputClusterNumber <- shiny::reactive({
        txt <- input$clusterNumbers
        if (is.null(txt) || !nzchar(txt)) {
            return(integer(0))
        }
        vals <- .asIntegerQuiet(stringr::str_split(txt, ",")[[1]])
        vals[!is.na(vals)]
    }) %>% shiny::debounce(1000)
    invisible(st)
}

#' Sample Ids of the Two t-Test Groups
#'
#' @param metaD \code{metadata(sce)} list with \code{experiment_info}.
#' @param groupsVar Grouping column of \code{experiment_info}.
#' @param group1,group2 Levels of \code{groupsVar} forming the groups.
#'
#' @return A list with \code{group1} and \code{group2} sample ids, or
#'   \code{NULL} when \code{groupsVar} is not a valid column.
#'
#' @keywords internal
.groupsFromInput <- function(metaD, groupsVar, group1, group2) {
    info <- metaD$experiment_info
    if (is.null(groupsVar) || length(groupsVar) == 0 ||
        !groupsVar %in% colnames(info)) {
        return(NULL)
    }
    list(
        group1 = info[info[, groupsVar] %in% group1, "sample_id"],
        group2 = info[info[, groupsVar] %in% group2, "sample_id"]
    )
}

#' Build the Group Plot Reactives
#'
#' Adds the violin plots, the UpSet plot and the FlowSOM marker pie plot to
#' \code{st}.
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#'
#' @return Called for its side effect on \code{st}.
#'
#' @keywords internal
.buildGroupPlotReactives <- function(st, input) {
    st$violinPlotSelection <- shiny::reactive(input$violinSelection) %>%
        shiny::debounce(1000)
    st$vlnPlot <- shiny::reactive({
        input$applyName
        input$rmGrp
        outputList <- st$rv$outputList
        shiny::req(outputList)
        upsetSel <- input$upsetSelection
        if (length(upsetSel) < 3) upsetSel <- names(outputList)
        data <- .violinData(st$metaD[[st$somCodesName]], upsetSel,
            outputList, st$violinPlotSelection())
        .violinGgplot(data, "No groups with >= 2 SOM nodes")
    })
    st$VlnPlot2 <- shiny::reactive({
        input$applyName
        input$rmGrp
        outputList <- st$rv$outputList
        shiny::req(outputList)
        plotViolin2Func(st$sce, st$somCodesName, st$violinPlotSelection(),
            input$upsetSelection, outputList)
    })
    st$upSetPlot <- shiny::reactive({
        st$selectedUpdate()
        input$applyName
        input$rmGrp
        upsetPlotFunc(input$upsetSelection, st$rv$outputList, st$sce)
    })
    st$flowSOMPiePlot <- shiny::reactive({
        rs <- st$rsUsed()
        shiny::req(rs)
        .buildFlowSOMPiePlot(st$metaD[[st$somCodesName]], rs, st$colsUsed,
            ncol = input$flowSOMPieCols, maxPies = input$flowSOMPieMax)
    })
    invisible(st)
}

#' Build the Node-Properties Reactives
#'
#' Per-node ambiguity of the subsample (computed once per session and
#' cached) and the node-properties data frame joining \code{SOM_stats},
#' the known-label node summary (when present) and the ambiguity summary.
#'
#' @inheritParams .buildDimRedPlotReactives
#'
#' @return Called for its side effect of adding reactives to \code{st}.
#'
#' @keywords internal
.buildNodePropsReactives <- function(st, input) {
    st$nodeAmbiguity <- shiny::reactive({
        nodeAmbiguity(st$sce_subsampled, somCodesName = st$somCodesName,
            markers = st$colsUsed)
    }) %>% shiny::bindCache(st$sce_subsampled)
    st$nodePropsDf <- shiny::reactive({
        stats <- S4Vectors::metadata(st$sce)$SOM_stats
        shiny::req(stats)
        labels <- st$metaD$node_labels
        df <- stats
        if (!is.null(labels)) {
            df <- dplyr::left_join(df, labels, by = "id")
        }
        # ambiguity also has an 'n' (cells sampled); keep the SOM_stats one
        amb <- st$nodeAmbiguity()
        if (!is.null(amb)) {
            df <- dplyr::left_join(df, amb, by = "id", suffix = c("", ".amb"))
        }
        df[order(df$id), , drop = FALSE]
    }) %>% shiny::bindCache(st$sce_subsampled)
    invisible(st)
}
