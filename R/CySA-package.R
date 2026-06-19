#' CySA: Interactive Cluster Selector for Cytometry Data
#'
#' Provides an interactive Shiny application for selecting and visualizing
#' clusters from flow cytometry data stored in \code{SingleCellExperiment}
#' objects.
#'
#' @keywords internal
#' @import shiny
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar sidebarMenu menuItem dashboardBody
#' @import shinyjqui
#' @import ggplot2
#' @import dplyr
#' @import stringr
#' @import ComplexHeatmap
#' @import dendextend
#' @import viridis
#' @import cowplot
#' @import ggplotify
#' @import grid
#' @importFrom shinyjs useShinyjs extendShinyjs
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom data.table rbindlist
#' @importFrom purrr map2 set_names imap is_empty
#' @importFrom plotly renderPlotly plotlyOutput ggplotly layout event_register event_data add_trace config
#' @importFrom DT renderDT DTOutput
#' @importFrom collapsibleTree collapsibleTreeOutput renderCollapsibleTree
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom S4Vectors metadata
#' @importFrom CATALYST channels cluster_codes cluster_ids
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap umap.defaults
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette pdf dev.off recordPlot replayPlot
#' @importFrom stats prcomp t.test sd reorder
#' @importFrom utils globalVariables str
#' @importFrom Matrix rowSums
#' @importFrom graphics points
#' @importFrom raster plot
CySA_default_cluster_cols <- function() {
  c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
    "#666666", "#999999"
  )
}

"_PACKAGE"
