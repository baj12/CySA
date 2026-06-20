
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CySA

<!-- badges: start -->
<!-- badges: end -->

**CySA** provides an interactive Shiny application for selecting and visualizing
clusters from flow-cytometry data stored in
[`SingleCellExperiment`](https://bioconductor.org/packages/SingleCellExperiment)
objects. It is designed to work with SOM-based clustering outputs such as those
produced by [FlowSOM](https://bioconductor.org/packages/FlowSOM) and curated by
the [CATALYST](https://bioconductor.org/packages/CATALYST) workflow.

## Installation

Install the development version from GitHub with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# After Bioconductor acceptance:
BiocManager::install("CySA")

# Development version:
BiocManager::install("bernd/CySA")  # replace with your GitHub user/repo
```

## Example

``` r
library(CySA)
library(SingleCellExperiment)

# sce should contain SOM_codes and SOM_stats in metadata(sce)
prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 10000)

# Launch the interactive selector
cs <- clusterSelector(
  sce = prepped$sce,
  sce_subsampled = prepped$sce_subsampled,
  dList = prepped$dList
)

shiny::shinyApp(ui = cs[[1]], server = cs[[2]])
```

## Citation

If you use CySA in your research, please cite the package and the relevant
Bioconductor workflows.
