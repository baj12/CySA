# ClusterSelector Shiny App - Package Extraction Analysis

**Date:** 2026-06-17  
**Author:** Claude Code  
**Purpose:** Analyze what needs to be transferred to extract the ClusterSelector Shiny app into a standalone R package

---

## 1. Problem Statement

The `clusterSelector.shiny.R` file contains a complex Shiny application for interactive cluster selection and visualization in flow cytometry data analysis. Currently it:
- Lives in the CyDa package root directory (not in `R/`)
- Sources external function files (`shinyAppFunctions.R`)
- Has extensive dependencies on CyDa internal functions
- Uses hardcoded absolute paths

**Goal:** Extract this into a standalone R package that can be installed and used independently.

---

## 2. Current Architecture

### 2.1 Main Files

| File | Location | Purpose |
|------|----------|---------|
| `clustereSelector.shiny.R` | Package root | Main Shiny app (1907 lines) |
| `prep.clusterSelector.shiny.R` | Package root | Pre-computation script for dimSelection |
| `shinyAppFunctions.R` | Package root | Helper functions for plots and stats |

### 2.2 How the App is Called

Example from `dendrogramApp.3.paloma.R`:

```r
# 1. Load SCE object with SOM codes
sce <- read_rds(somFile)

# 2. Create subsampled version for performance
source("prep.clusterSelector.shiny.R")  # Creates sce_subsampled

# 3. Compute dendrogram from SOM codes
hc = hclust(dist(somCodes))
dend <- as.dendrogram(hc)

# 4. Build dendTable for tree navigation
dendTable = data.frame("parent"=NA,"child"="rt", ...)
recDendTable(dend)

# 5. Create cluster-patient contingency table
clusterPatientTable = table(colData(sce)[,c("sample_id", "cluster_id")])

# 6. Build SOM raster for visualization
somRasterData <- data.frame(x = ..., y = ..., somCodes)
somRasterObj <- rasterFromXYZ(xyz = somRasterData)

# 7. Define marker pairs for 2D plots
dList = list(
  d1 = c("FSC.A", "SSC.A"),
  d2 = c("FSC.A", "FSC.H"),
  ...
)

# 8. Launch app
cS = clusterSelector(sce, sce_subsampled, outputList, dList, dend, dendTable, 
                     clusterPatientTable, somRasterData, somRasterObj, ...)
shinyApp(ui = cS[[1]], server = cS[[2]], options = list(launch.browser = T))
```

---

## 3. Input Specification

### 3.1 Function Signature

```r
clusterSelector <- function(
  sce,                    # SingleCellExperiment - MAIN input
  sce_subsampled,         # SingleCellExperiment - subsampled for performance
  outputList = list(),    # list of named cluster groupings
  colTree = NULL,         # collapsible tree object (optional)
  dList,                  # list of 2D plot axis pairs
  dend,                   # dendrogram object
  dendTable,              # data.frame for dendrogram navigation
  clusterPatientTable,    # table: sample_id x cluster_id
  somCodesName = "SOM_codes",  # metadata slot name
  nPlots = 6,             # number of 2D SOM plots
  somRasterData,          # data.frame for SOM raster
  somRasterObj,           # raster object
  env = environment()
)
```

### 3.2 Required SCE Metadata

The `SingleCellExperiment` object must contain:

**In `metadata(sce)`:**
- `SOM_codes` (or custom `somCodesName`): Matrix of SOM cluster codes
- `SOM_stats`: Data frame with columns:
  - `median`, `mean`, `rdQu`, `max`, `n`, `id`
- `experiment_info`: Data frame with sample metadata including:
  - `sample_id` (required)
  - `group` (optional, for grouping)
  - Numeric columns for stats comparison
  - Factor columns for t-test grouping

**In `colData(sce)`:**
- `cluster_id`: Factor assigning cells to SOM clusters
- `sample_id`: Factor assigning cells to samples

**In `assays(sce)`:**
- At least one assay (typically "exprs") with markers as rows, cells as columns

### 3.3 Input Object Details

#### `outputList`
```r
list(
  "ClusterName1" = c(1, 5, 12, 23),    # SOM cluster IDs
  "ClusterName2" = c(3, 7, 8),
  "Rest" = c(...)                       # automatically computed
)
```

#### `dList`
```r
list(
  d1 = c("FSC.A", "SSC.A"),    # marker names for plot 1 (X, Y)
  d2 = c("FSC.A", "FSC.H"),    # marker names for plot 2
  ...
)
```

#### `dendTable`
```r
data.frame(
  parent = c(NA, "rt", "rt", "rt.1", ...),
  child = c("rt", "rt.1", "rt.2", "rt.1.1", ...),
  nleaf = c(40, 20, 20, 10, ...),
  indexString = c("", "[[1]]", "[[2]]", "[[1]][[1]]", ...)
)
```

#### `somRasterData`
```r
data.frame(
  x = c(1,2,3,4,5,6,7,8,...),
  y = c(1,1,1,1,2,2,2,2,...),
  CD31.A = c(...),    # SOM code values for each grid position
  Dead.A = c(...),
  ...
)
```

---

## 4. Dependencies Analysis

### 4.1 R Package Dependencies (NAMESPACE)

```r
# Core
Imports:
  shiny,
  shinydashboard,
  shinydashboardPlus,
  shinyjs,
  plotly,
  ggplot2,
  ggthemr,
  ggplotify,
  jqr (for jqui_resizable),
  sortable,
  
# Data handling
  SingleCellExperiment,
  S4Vectors,
  flowCore,
  SummarizedExperiment,
  dplyr,
  tidyr,
  tibble,
  data.table,
  stringr,
  readr,
  purrr,
  
# Visualization
  ComplexHeatmap,
  dendextend,
  viridis,
  RColorBrewer,
  grid,
  gridExtra,
  
# Dimensionality reduction
  Rtsne,
  umap,
  
# Utilities
  DT,
  reshape2,
  rlang
```

### 4.2 CyDa Functions Used

These functions are currently called from CyDa's R/ directory:

| Function | Source File | Purpose |
|----------|-------------|---------|
| `plotSOMScatter()` | `R/plotSOMScatter.R` | Plot SOM grid with 2 markers |
| `plotScatterBJ()` | `R/plotScatterBJ.R` | Scatter plot for flow cytometry |
| `channels()` | `R/channels.R` | Get channel names from SCE |

### 4.3 Functions from `shinyAppFunctions.R`

These 11 functions must be transferred:

| Function | Purpose |
|----------|---------|
| `highlight_df()` | Extract SOM code coordinates for highlighting |
| `countBarPlotFunc()` | Create count bar plot data |
| `PercentBarPlotFunc()` | Create percentage bar plot data |
| `ggsomPlot()` | Add highlighting to SOM plot |
| `somPlot()` | Render interactive SOM scatter plot |
| `tsneFunc()` | Compute t-SNE embedding |
| `plotViolinFunc()` | Create violin plot data (marker view) |
| `plotViolin2Func()` | Create violin plot data (group view) |
| `drawProjection()` | Draw 2D projection with highlighting |
| `upsetPlotFunc()` | Create UpSet plot for set overlaps |
| `rasterFromXYZ()` | Convert to raster (currently from `raster` pkg) |

---

## 5. UI Components

### 5.1 Sidebar (Parameters)

| Input ID | Type | Purpose |
|----------|------|---------|
| `tabs` | tabMenu | Navigation |
| `d2Axes` | selectizeInput | Choose plot index |
| `d{1-6}.1`, `d{1-6}.2` | selectInput | X/Y axis markers per plot |
| `selectMode` | radioButtons | view/add/remove selection mode |
| `samples2plot` | selectizeInput | Samples to display |
| `clusterNumbers` | textInput | Manual cluster ID entry |
| `clusterName` | textInput | Name for selection |
| `clusterNameSelect` | selectizeInput | Select named groups |
| `clusterNameRM` | selectizeInput | Group to remove |
| `groupRM` | selectizeInput | Groups to remove (multi) |
| `colorbyGroups` | selectizeInput | Groups for coloring |
| `dimRedSelection` | selectizeInput | Markers for dim reduction |
| `perplexity` | numericInput | t-SNE perplexity |
| `n_neighbors` | numericInput | UMAP neighbors |
| `compareStatsTo` | selectizeInput | Stats comparison variable |
| `relativeTo` | selectizeInput | Parent population for % |
| `singleNode` | numericInput | Single SOM node counts |
| `groupsVar` | selectizeInput | Factor for t-test |
| `group1`, `group2` | selectizeInput | Groups for t-test |

### 5.2 Main Body Outputs

| Output ID | Type | Content |
|-----------|------|---------|
| `plot` | collapsibleTree | Interactive dendrogram |
| `scatter` | plotly | 2D scatter plot |
| `dend` | plot | Dendrogram with highlighting |
| `somData{1-6}` | plotly | SOM 2D projections |
| `tsne`, `umap`, `pca` | plotly | Dimensionality reduction |
| `somRaster` | plot | SOM raster overview |
| `somRasterSelect` | plotly | Raster cluster selection |
| `VlnPlot`, `VlnPlot2` | plot | Violin plots |
| `UpSet` | plot | UpSet plot |
| `cellCounts` | DT | Cell count table |
| `cellPercentages` | renderPrint | Percentage stats |
| `CountBar`, `PercentBar` | plot | Bar visualizations |
| `ttestResult` | renderPrint | T-test output |
| `somClusters` | renderText | Selected cluster IDs |

---

## 6. Reactive Logic

### 6.1 Key Reactives

| Reactive | Type | Purpose |
|----------|------|---------|
| `rsUsed` | reactiveVal | Currently selected cluster IDs |
| `dimSelection` | reactiveVal | Axis selections and zoom state |
| `selectedPoints` | reactiveVal | Points from plotly selection |
| `activePlot` | reactiveVal | Which 2D plot is active |
| `groupsInput` | reactive | Debounced group selection for t-test |
| `tsne`, `umap`, `pca` | reactive | Dimensionality reduction (debounced) |
| `violinPlotSelection` | reactive | Selected markers for violin |
| `dimRedSelection` | reactive | Markers for dim reduction |

### 6.2 Event Observers

| Event | Action |
|-------|--------|
| `plotly_selected` (somData1-6) | Update selected clusters |
| `plotly_selected` (tsne/umap/pca) | Update selected clusters |
| `plotly_selected` (scatterPlot) | Box selection on scatter |
| `plotly_selected` (somGrid) | Raster cluster selection |
| `plotly_relayout` | Handle zoom |
| `input$node` | Dendrogram node click |
| `input$applyName` | Save named selection |
| `input$rmGrp` | Remove group |
| `input$applyclusterNumbers` | Select by cluster ID |

---

## 7. Package Structure Proposal

```
ClusterSelector/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── clusterSelector.R          # Main function (extracted from .shiny.R)
│   ├── shinyAppFunctions.R        # Helper functions
│   ├── prepClusterSelector.R      # Preparation helper
│   ├── plotSOMScatter.R           # (copy from CyDa)
│   ├── plotScatterBJ.R            # (copy from CyDa)
│   └── utils.R                    # Shared utilities
├── man/                           # Roxygen documentation
├── inst/
│   └── app/                       # Shiny app structure (optional)
│       ├── ui.R
│       └── server.R
├── tests/
│   └── testthat/
│       └── test-clusterSelector.R
└── vignettes/
    └── clusterSelector.Rmd
```

---

## 8. Migration Steps

### Step 1: Create Package Skeleton
```r
usethis::create_package("ClusterSelector")
usethis::use_roxygen_md()
usethis::use_testthat()
```

### Step 2: Copy Core Files
- `clustereSelector.shiny.R` → `R/clusterSelector.R`
- `shinyAppFunctions.R` → `R/shinyAppFunctions.R`
- `prep.clusterSelector.shiny.R` → `R/prepClusterSelector.R`

### Step 3: Copy Dependencies from CyDa
- `plotSOMScatter.R`
- `plotScatterBJ.R`
- Any utility functions they depend on

### Step 4: Update Imports
- Replace `source()` calls with proper `@importFrom` roxygen tags
- Fix hardcoded paths (line 5 in clustereSelector.shiny.R)
- Remove `require()` calls, use `::` or imports

### Step 5: Fix the `source()` Line
Line 5 currently:
```r
source("/pasteur/helix/projects/scBiomarkers/bernd/cytometry/CyDa/shinyAppFunctions.R")
```
Should become:
```r
# Functions are now part of the package namespace
```

### Step 6: Export Function
Add to `NAMESPACE`:
```
export(clusterSelector)
export(prepareClusterSelectorData)
```

### Step 7: Documentation
Write roxygen2 documentation for:
- `clusterSelector()` - main function
- Input requirements
- Examples

### Step 8: Test
```r
devtools::load_all()
devtools::test()
```

---

## 9. Special Considerations

### 9.1 Hardcoded Paths
Line 5 in `clustereSelector.shiny.R`:
```r
source("/pasteur/helix/projects/scBiomarkers/bernd/cytometry/CyDa/shinyAppFunctions.R")
```
**Fix:** Remove - functions will be in package namespace.

### 9.2 Global Variables
The app uses `<<-` and `env` parameter for state:
```r
assign(x = "outputList", value = outputList, envir = env)
outputList = get("outputList", envir = env)
```
**Consider:** Use `shiny::reactiveValues` or return modified objects.

### 9.3 The `sce` Global
Many functions reference `sce` directly without it being passed:
```r
drawProjection <- function(df, rs, colorbyGroups, env=env){
  df$N = S4Vectors::metadata(sce)$SOM_stats$n  # sce not in scope!
```
**Fix:** Pass `sce` explicitly or use closure.

### 9.4 `rasterFromXYZ` Dependency
Currently uses `raster::rasterFromXYZ()`. Consider:
- Keep as dependency
- Or implement simple replacement (it's a small function)

### 9.5 `ggthemr` for Theming
```r
ggthemr('flat')
ggthemr("dust")
```
**Note:** `ggthemr` is archived on CRAN. Consider:
- Use `ggplot2::theme_minimal()` or similar
- Or suggest theme as parameter

---

## 10. Recommended Package Dependencies

```yaml
Imports:
  - R (>= 4.0.0)
  - shiny (>= 1.7.0)
  - shinydashboard
  - shinydashboardPlus
  - shinyjs
  - plotly
  - ggplot2
  - ggthemr (or suggest alternative)
  - ggplotify
  - jqr
  - sortable
  - SingleCellExperiment
  - S4Vectors
  - SummarizedExperiment
  - flowCore
  - dplyr
  - tidyr
  - tibble
  - data.table
  - stringr
  - purrr
  - ComplexHeatmap
  - dendextend
  - viridis
  - RColorBrewer
  - Rtsne
  - umap
  - DT
  - reshape2
  - raster (or implement rasterFromXYZ)

Suggests:
  - CyDa (for example data)
  - testthat
  - knitr
  - rmarkdown
```

---

## 11. Testing Strategy

### 11.1 Unit Tests
```r
test_that("highlight_df returns correct structure", {
  result <- highlight_df("CD3.A", "Dead.A", c(1,2,3), "SOM_codes", metaD)
  expect_equal(ncol(result), 3)
  expect_true(all(c("x", "y", "id") %in% names(result)))
})
```

### 11.2 Integration Test
```r
test_that("clusterSelector launches", {
  # Create minimal mock SCE
  sce <- mockSCE()
  result <- clusterSelector(sce, ...)
  expect_type(result, "list")
  expect_length(result, 2)  # ui, server
})
```

### 11.3 Visual Tests
- Compare plot outputs against known good references
- Test with example dataset from `CyDa`

---

## 12. Exported Functions

| Function | Export | Purpose |
|----------|--------|---------|
| `clusterSelector()` | Yes | Main Shiny app launcher |
| `prepareClusterSelectorData()` | Yes | Helper to prepare all inputs |
| `plotSOMScatter()` | Maybe | General-purpose plot function |
| `plotScatterBJ()` | Maybe | General-purpose plot function |

---

## 13. Example Usage (After Packaging)

```r
library(ClusterSelector)

# Load your SCE object
sce <- readRDS("path/to/sce_with_som.RDS")

# Prepare data (helper function)
prepared <- prepareClusterSelectorData(
  sce = sce,
  n_subsample = 100000,
  dList = list(
    d1 = c("FSC.A", "SSC.A"),
    d2 = c("CD3.A", "CD19.A")
  )
)

# Launch app
shinyApp(
  ui = prepared$ui,
  server = prepared$server
)

# Or simpler:
runClusterSelector(sce)  # wrapper that does all preparation
```

---

## 14. Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Breaking changes in CyDa functions | Copy functions, don't depend on CyDa |
| `ggthemr` removed from CRAN | Use alternative theming, make configurable |
| Large memory footprint | Document requirements, add subsampling |
| Complex input requirements | Provide `prepareClusterSelectorData()` helper |
| Hardcoded paths | Audit all paths, use `system.file()` |

---

## 15. Recommended Next Steps

1. **Create package skeleton** with `usethis::create_package()`
2. **Copy and adapt** the three main source files
3. **Copy dependencies** (`plotSOMScatter.R`, `plotScatterBJ.R`)
4. **Write NAMESPACE** with proper imports
5. **Create helper function** `prepareClusterSelectorData()` that:
   - Creates subsampled SCE
   - Computes dendrogram and dendTable
   - Builds somRasterData/Obj
   - Creates clusterPatientTable
   - Sets up default dList
6. **Write documentation** with examples
7. **Test with real data** from CyDa examples
8. **Consider renaming** to avoid confusion (e.g., `FlowClusterSelector`)

---

## 16. File Size Estimates

| Component | Lines | Priority |
|-----------|-------|----------|
| clusterSelector (main) | ~1900 | Core |
| shinyAppFunctions | ~410 | Core |
| prep.clusterSelector | ~100 | Helper |
| plotSOMScatter | ~150 | Dependency |
| plotScatterBJ | ~200 | Dependency |
| **Total** | **~2760** | |

---

## 17. Conclusion

The ClusterSelector Shiny app is a **self-contained but complex** application that can be extracted into a standalone package. Key observations:

1. **Well-structured inputs**: The function signature clearly defines what's needed
2. **Externalized helpers**: `shinyAppFunctions.R` is already separated
3. **Minimal CyDa dependencies**: Only 2-3 plotting functions are needed
4. **State management**: Uses environment for `outputList` - could be improved
5. **Hardcoded paths**: One critical path to fix (line 5)

**Estimated effort:** 2-3 days for initial package, 1 week for polished release with tests and vignettes.

---

## 18. Tools to Create the New Package from Scratch

Several R packages help create new R packages from scratch. The most commonly used are:

### 18.1 `usethis` (Recommended)

The modern standard for package development. Useful functions for the ClusterSelector package:

```r
# Create a new package directory
usethis::create_package("../CySA")

# Navigate into the package (interactive)
usethis::proj_activate("path/to/ClusterSelector")

# Set up infrastructure
usethis::use_roxygen_md()          # Use roxygen2 with markdown
usethis::use_testthat()            # Add testthat infrastructure
usethis::use_mit_license()         # Add a license
usethis::use_readme_rmd()          # Create README.Rmd
usethis::use_vignette("intro")     # Create a vignette
usethis::use_pkgdown()             # Website (optional)

# Add dependencies
usethis::use_package("shiny")
usethis::use_package("shinydashboard")
usethis::use_package("plotly")
# ... etc

# Create R files
usethis::use_r("clusterSelector")
usethis::use_r("shinyAppFunctions")
```

### 18.2 `devtools`

Higher-level wrapper around many package tasks. Typically used together with `usethis`:

```r
devtools::create("path/to/ClusterSelector")    # Alternative to usethis
devtools::document()                           # Generate NAMESPACE/man files
devtools::check()                              # Run R CMD check
devtools::install()                            # Install locally
devtools::test()                               # Run tests
devtools::build()                              # Build source tarball
devtools::load_all()                           # Load package for testing
```

### 18.3 `pkgdown` (Optional)

For creating documentation websites:

```r
usethis::use_pkgdown()
pkgdown::build_site()
```

### 18.4 `available`

Check if a package name is available on CRAN/Bioconductor:

```r
available::available("ClusterSelector")
```

### 18.5 `biocthis`

If you want to create a Bioconductor-style package:

```r
biocthis::use_bioc_pkg_templates()
```

### 18.6 Recommended Workflow for ClusterSelector

```r
# Step 1: Create package
usethis::create_package("~/ClusterSelector")

# Step 2: Set up infrastructure
usethis::use_roxygen_md()
usethis::use_testthat()
usethis::use_mit_license("Your Name")
usethis::use_readme_rmd()

# Step 3: Add dependencies
usethis::use_package("shiny")
usethis::use_package("shinydashboard")
usethis::use_package("shinydashboardPlus")
usethis::use_package("shinyjs")
usethis::use_package("plotly")
usethis::use_package("ggplot2")
usethis::use_package("SingleCellExperiment")
usethis::use_package("S4Vectors")
usethis::use_package("dplyr")
usethis::use_package("tidyr")
usethis::use_package("ComplexHeatmap")
usethis::use_package("Rtsne")
usethis::use_package("umap")
# ... add remaining dependencies

# Step 4: Copy R files into R/ directory
# Copy: clustereSelector.shiny.R, shinyAppFunctions.R,
#       prep.clusterSelector.shiny.R,
#       plotSOMScatter.R, plotScatterBJ.R

# Step 5: Generate docs and check
devtools::document()
devtools::check()
devtools::install()
devtools::test()
```

### 18.7 Notes on Dependency Versioning

When using `usethis::use_package()`:
- It adds packages to `DESCRIPTION` under `Imports:` by default
- It does NOT add them to `NAMESPACE`
- You must use `@importFrom pkg function` in roxygen comments or `pkg::function()` in code
- For the ClusterSelector package, prefer `pkg::function()` style to keep NAMESPACE minimal

### 18.8 Conclusion on Tool Choice

For this extraction task, the best combination is:

- **Scaffolding:** `usethis::create_package()`
- **Iterative development:** `devtools::document()`, `devtools::check()`, `devtools::test()`
- **Documentation website:** `pkgdown` (optional, after package is functional)
- **Name checking:** `available::available()` before committing to a package name

---

*End of analysis document.*
