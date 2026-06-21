# CySA 0.99.3

## Bioconductor submission fixes

* Wrapped the second occurrence of a long `outputList[["Rest"]]` assignment.

# CySA 0.99.2

## Bioconductor submission fixes

* Wrapped remaining long lines flagged by `BiocCheck`.

# CySA 0.99.1

## Bioconductor submission fixes

* Replaced `<<-` with `shiny::reactiveVal()` for app state.
* Removed debug `cat()`/`print()` calls; kept `print()` only inside PDF device
  context.
* Applied `styler` formatting and wrapped selected long lines.

# CySA 0.99.0

## New features

* Initial Bioconductor submission of `CySA`, an interactive Shiny dashboard for
  selecting and visualizing clusters from flow-cytometry data stored in
  `SingleCellExperiment` objects.
* `clusterSelector()` returns a Shiny app object that can be launched with
  `shiny::runApp()`.
* `prepClusterSelectorData()` subsamples a `SingleCellExperiment` and builds the
  inputs required by `clusterSelector()`.
* `plotSOMScatter()` and `plotScatterBJ()` provide ggplot2-based SOM and scatter
  visualizations.

## Bioconductor readiness

* Replaced broad package-level `@import` directives with targeted
  `@importFrom` roxygen tags to reduce namespace masking.
* Removed the archived `ggthemr` dependency; theming now uses `theme_minimal()`
  and explicit `scale_*_manual()` scales.
* Added `NEWS.md`, expanded `README.Rmd`, and a package vignette.
