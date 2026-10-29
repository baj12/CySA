# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

utils::globalVariables(c(
    ".", ".data", "..ncount..", "counts", "groups", "id", "marker", "expr",
    "grpName", "Percent", "n", "value", "variable", "cluster", "colGrp",
    "x", "y", "xend", "yend", "customdata", "label", "tsne1", "tsne2",
    "umap1", "umap2", "xe", "ye",
    "pc1", "pc2", "rowElement", "cellCounts", "cellPercentages",
    "dend", "dendTable", "clusterPatientTable", "somRasterData",
    "somRasterObj", "sce_subsampled", "metaD", "sceRN", "sceCN",
    # Names loaded from the saved prep bundle by .applyLabelStatsToBundle().
    "dList", "dend_obj", "dend_table",
    "chx", "nNsub", "sce_subsampledRN", "df", "dfPlot",
    "outputList", "dimSelection", "activePlot", "rsUsed",
    "inputClusterNumber", "violinPlotSelection", "groupsInput",
    "countBarPlot", "PercentBarPlot", "dendPlot", "tsnePlot",
    "umapPlot", "pcaPlot", "scatterPlot", "somRasterPlot", "flowSOMPiePlot",
    "vlnPlot", "VlnPlot2", "upSetPlot", "selectedUpdate", "selectedUpdate2",
    "sample2PlotDb", "choicesRV",
    "sample_id", "group", "val", "somNode", "N", "thrdQu",
    "scatterPercentile", "somColorVar", "somSizeVar",
    ".dens_col",
    "rsUsed_d",
    ".pt", ".stroke",
    # aes() columns of the marker ridge trace plot (plot-ridge-trace.R).
    "height", "fillKey", "lineColour", "datatype", "ridgeGroup",
    "linewidth", "cellId", "isExtreme", "point_shape", "point_size",
    "point_stroke",
    # aes() columns and testServer bare names of the diagnostics, ordered
    # cluster heatmap and diffcyt panels (plot-diagnostics.R,
    # diffcyt-helpers.R, server-reactives-diagnostics.R).
    "mds1", "mds2", "nrs", "color_by", "cluster_id",
    "diagMdsPlot", "diagExprHeatmapPlot", "diagNrsPlot",
    "orderedHeatmapPlot", "diffResult",
    # aes() column and output name of the cluster-ambiguity diagnostic.
    "pctAmbiguous", "diagAmbiguityPlot",
    # aes() columns of the node-properties 2D plot (node-properties.R,
    # server-reactives.R) and its output name.
    "purity", "cell_type_majority", "node_dispersion",
    "meanDistToCode", "medianRatio", "nodePropsDf", "nodeAmbiguity"
))

#' Cluster Selector Shiny Application
#'
#' Creates an interactive Shiny dashboard for selecting and visualizing
#' clusters from a SingleCellExperiment object.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing the full dataset.
#' @param sce_subsampled A subsampled \code{SingleCellExperiment} for
#'   performance-sensitive plots.
#' @param outputList A named list of cluster groupings.
#' @param colTree Optional collapsible tree object.
#' @param dList List of marker pairs for 2D plots.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param nPlots Number of 2D SOM plots to display.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somRasterObj Raster object for SOM visualization. If \code{NULL},
#'   the SOM grid dimensions are inferred from \code{somRasterData}.
#' @param fsom Optional \code{FlowSOM} object. If provided, an interactive
#'   FlowSOM star plot is shown.
#' @param env Environment used to store mutable state (legacy argument).
#' @param verbose Logical indicating whether to show detailed runtime messages.
#'
#' @return A \code{\link[shiny]{shinyApp}} object. Launch the app with
#'   \code{shiny::runApp(app)}.
#'
#' @examples
#' sce <- CySA_example_sce()
#' prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
#' som_codes <- S4Vectors::metadata(sce)$SOM_codes
#' dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
#' dendTable <- data.frame(
#'     id = seq_len(nrow(som_codes)),
#'     label = rownames(som_codes),
#'     stringsAsFactors = FALSE
#' )
#' clusterPatientTable <- table(
#'     sample_id = sce$sample_id,
#'     cluster_id = sce$cluster_id
#' )
#' markers <- S4Vectors::metadata(sce)$map$colsUsed
#' somRasterData <- data.frame(
#'     x = rep(seq_len(5), length.out = nrow(som_codes)),
#'     y = rep(seq_len(2), each = nrow(som_codes) / 2),
#'     id = seq_len(nrow(som_codes))
#' )
#' for (m in markers) {
#'     somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
#' }
#' arr <- array(
#'     data = seq_len(10 * 10 * length(markers)),
#'     dim = c(10, 10, length(markers))
#' )
#' somRasterObj <- raster::brick(arr)
#' names(somRasterObj) <- markers
#'
#' app <- clusterSelector(
#'     sce = prepped$sce,
#'     sce_subsampled = prepped$sce_subsampled,
#'     dList = prepped$dList,
#'     dend = dend,
#'     dendTable = dendTable,
#'     clusterPatientTable = clusterPatientTable,
#'     somRasterData = somRasterData,
#'     somRasterObj = somRasterObj
#' )
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
#'
#' @export
clusterSelector <- function(sce,
                            sce_subsampled,
                            outputList = list(),
                            colTree = NULL,
                            dList,
                            dend,
                            dendTable,
                            clusterPatientTable,
                            somCodesName = "SOM_codes",
                            nPlots = 6,
                            somRasterData,
                            somRasterObj,
                            fsom = NULL,
                            env = environment(),
                            verbose = FALSE) {
    csa <- .buildClusterSelectorApp(
        sce = sce, sce_subsampled = sce_subsampled, outputList = outputList,
        colTree = colTree, dList = dList, dend = dend,
        dendTable = dendTable, clusterPatientTable = clusterPatientTable,
        somCodesName = somCodesName, nPlots = nPlots,
        somRasterData = somRasterData, somRasterObj = somRasterObj,
        fsom = fsom, env = env, verbose = verbose
    )
    shiny::shinyApp(ui = csa$ui, server = csa$server)
}
