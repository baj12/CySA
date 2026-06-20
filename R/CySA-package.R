# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding assistant.
# All code is redistributed under the package LICENSE.

#' CySA: Interactive Cluster Selector for Cytometry Data
#'
#' Provides an interactive Shiny application for selecting and visualizing
#' clusters from flow cytometry data stored in \code{SingleCellExperiment}
#' objects.
#'
#' @keywords internal
#'
#' @import shiny
#' @import ggplot2
#'
#' @importFrom shinydashboard dashboardBody
#' @importFrom shinydashboard dashboardHeader
#' @importFrom shinydashboard dashboardPage
#' @importFrom shinydashboard dashboardSidebar
#' @importFrom shinydashboard menuItem
#' @importFrom shinydashboard sidebarMenu
#' @importFrom shinyjs useShinyjs extendShinyjs
#'
#' @importFrom dplyr `%>%`
#' @importFrom dplyr as_tibble
#' @importFrom dplyr count
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr join_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr pull
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom stringr str_replace
#' @importFrom stringr str_split
#' @importFrom tidyr all_of
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom data.table rbindlist
#' @importFrom purrr map2
#' @importFrom purrr set_names
#' @importFrom purrr imap
#' @importFrom purrr is_empty
#'
#' @importFrom plotly add_trace
#' @importFrom plotly config
#' @importFrom plotly event_data
#' @importFrom plotly event_register
#' @importFrom plotly ggplotly
#' @importFrom plotly layout
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom DT DTOutput
#' @importFrom DT renderDT
#' @importFrom collapsibleTree collapsibleTreeOutput
#' @importFrom collapsibleTree renderCollapsibleTree
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom grid unit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' @importFrom viridis scale_fill_viridis
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom SingleCellExperiment int_colData
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom CATALYST channels
#' @importFrom CATALYST cluster_codes
#' @importFrom CATALYST cluster_ids
#'
#' @importFrom Matrix rowSums
#' @importFrom raster plot
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom umap umap.defaults
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics points
#' @importFrom stats prcomp
#' @importFrom stats quantile
#' @importFrom stats reorder
#' @importFrom stats sd
#' @importFrom stats t.test
#' @importFrom utils globalVariables
#' @importFrom utils head
#' @importFrom utils str
"_PACKAGE"

#' Default cluster color palette
#'
#' Returns the default 20-color palette used by CySA for cluster
#' visualizations.
#'
#' @return A character vector of hex colors.
#'
#' @examples
#' CySA_default_cluster_cols()
#'
#' @export
CySA_default_cluster_cols <- function() {
  c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
    "#666666", "#999999"
  )
}
