# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# utils_clusterSelector.R ----
# Internal app builder. Validates inputs, pre-computes cached objects, and
# returns the UI/server pair for clusterSelector().

#' Build clusterSelector UI and Server
#'
#' Creates the UI/server pair for the clusterSelector Shiny app.
#'
#' @inheritParams clusterSelector
#'
#' @return A named list with elements \code{ui} and \code{server}.
#'
#' @keywords internal
#' @name INTERNAL_buildClusterSelectorApp
.buildClusterSelectorApp <- function(sce, sce_subsampled, outputList = list(),
    colTree = NULL, dList, dend, dendTable, clusterPatientTable,
    somCodesName = "SOM_codes", nPlots = 6, somRasterData, somRasterObj,
    fsom = NULL, env = environment(), verbose = FALSE) {
    # nocov start
    metaD <- .validateClusterSelectorInputs(
        sce, sce_subsampled, outputList, dList, dend, dendTable,
        clusterPatientTable, somRasterData, somCodesName, nPlots, fsom
    )
    if ("outputList" %in% ls(envir = env)) {
        outputList <- get("outputList", envir = env)
    }
    clusterLevels <- if (is.factor(sce$cluster_id)) {
        levels(sce$cluster_id)
    } else {
        sort(unique(as.integer(sce$cluster_id)))
    }
    outputList <- .initializeOutputList(outputList, clusterLevels)
    colsUsed <- metaD$map$colsUsed
    raster <- .subsetSomRaster(somRasterData, somRasterObj, colsUsed)
    assign(x = "outputList", value = outputList, envir = env)
    ui <- .buildClusterSelectorUI(sce, outputList, colTree, dList, metaD,
        somCodesName, nPlots, fsom)
    server <- .buildClusterSelectorServer(
        sce = sce, sce_subsampled = sce_subsampled,
        outputList = outputList, colTree = colTree, dList = dList,
        dend = dend, dendTable = dendTable,
        clusterPatientTable = clusterPatientTable,
        somCodesName = somCodesName, nPlots = nPlots,
        somRasterData = raster$data, somRasterObj = raster$obj,
        fsom = fsom, env = env, verbose = verbose
    )
    list(ui = ui, server = server)
    # nocov end
}

#' Restrict the SOM Raster Inputs to the Clustering Markers
#'
#' @param somRasterData Data frame with \code{x}, \code{y}, optional
#'   \code{id} and one column per marker.
#' @param somRasterObj Raster brick/stack, array, or \code{NULL}.
#' @param colsUsed Clustering markers.
#'
#' @return A list with the subset \code{data} (with an \code{id} column) and
#'   \code{obj}.
#'
#' @keywords internal
.subsetSomRaster <- function(somRasterData, somRasterObj, colsUsed) {
    if (!"id" %in% colnames(somRasterData)) {
        somRasterData$id <- seq_len(nrow(somRasterData))
    }
    keepCols <- intersect(c("x", "y", "id", colsUsed), colnames(somRasterData))
    if (inherits(somRasterObj, c("RasterBrick", "RasterStack"))) {
        somRasterObj <- raster::subset(somRasterObj, colsUsed)
    } else if (is.array(somRasterObj)) {
        layers <- which(dimnames(somRasterObj)[[3]] %in% colsUsed)
        somRasterObj <- somRasterObj[, , layers, drop = FALSE]
    }
    list(data = somRasterData[, keepCols], obj = somRasterObj)
}

#' Build the clusterSelector Dashboard UI
#'
#' @inheritParams clusterSelector
#' @param metaD \code{metadata(sce)} list.
#'
#' @return A \code{shinydashboard::dashboardPage}.
#'
#' @keywords internal
.buildClusterSelectorUI <- function(sce, outputList, colTree, dList, metaD,
    somCodesName, nPlots, fsom) {
    rnSCE <- rownames(sce)
    colsUsed <- metaD$map$colsUsed
    shiny::tagList(
        .buildFaviconLink(),
        shinydashboard::dashboardPage(
            .buildDashboardHeader(),
            sidebar = .buildClusterSelectorSidebar(sce, outputList, nPlots),
            body = shinydashboard::dashboardBody(
                .buildFirstBodyRow(colTree, nPlots),
                .buildSOM2DPlotsBox(nPlots, colsUsed, outputList,
                    markers = rnSCE, dList = dList),
                .buildSixStaticSOMBox(markers = rnSCE, dList = dList),
                .buildFlowSOMPieBox(),
                .buildFlowSOMStarsBox(fsom),
                .buildStatsBox(metaD, somCodesName, outputList),
                .buildSOMRasterBox(),
                .buildViolinBox(metaD[[somCodesName]], colsUsed),
                .buildRidgeBox(metaD[[somCodesName]], colsUsed),
                .buildUpSetBox(outputList),
                .buildDiagnosticsBox(metaD, outputList),
                .buildOrderedHeatmapBox(outputList,
                    colnames(metaD[[somCodesName]])),
                .buildDiffTestBox(metaD, outputList, rnSCE)
            )
        )
    )
}

#' Browser-Tab Favicon
#'
#' Four blue SOM nodes around one red "selected" node -- the same red/blue
#' convention used throughout the app (e.g. the interactive dendrogram and
#' 2D plots) -- as an inline SVG data URI, so the icon needs no extra
#' package resource file or \code{shiny::addResourcePath()} registration
#' and works the same way in any deployment.
#'
#' @return A \code{shiny::tags$head()} tag holding the favicon
#'   \code{<link>}.
#'
#' @keywords internal
.buildFaviconLink <- function() {
    svg <- paste0(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32">',
        '<circle cx="8" cy="8" r="6.5" fill="#357ca5"/>',
        '<circle cx="24" cy="8" r="6.5" fill="#357ca5"/>',
        '<circle cx="8" cy="24" r="6.5" fill="#357ca5"/>',
        '<circle cx="24" cy="24" r="6.5" fill="#357ca5"/>',
        '<circle cx="16" cy="16" r="7" fill="#dc3545"/>',
        '</svg>'
    )
    href <- paste0("data:image/svg+xml,",
        utils::URLencode(svg, reserved = TRUE))
    shiny::tags$head(
        shiny::tags$link(rel = "icon", type = "image/svg+xml", href = href)
    )
}

#' Dashboard Header With the Running CySA Version
#'
#' The version tag lets users and bug reports reference the right build.
#' The AdminLTE logo box (.main-header .logo) is 50px high with an
#' inherited 50px line-height and overflow:hidden, which would push a second
#' line out of view; the embedded style relaxes the line-height so both the
#' title and the version fit. The version sits on its own line
#' (display:block) because the fixed-width logo box clips inline suffixes.
#'
#' @return A \code{shinydashboard::dashboardHeader}.
#'
#' @keywords internal
.buildDashboardHeader <- function() {
    csaVersion <- tryCatch(
        as.character(utils::packageVersion("CySA")),
        error = function(e) "devel"
    )
    shinydashboard::dashboardHeader(
        title = shiny::span(
            shiny::tags$style(shiny::HTML(paste(
                ".main-header .logo { line-height: 1.2; padding-top: 4px;",
                "padding-bottom: 4px; }"
            ))),
            "Cluster Selector",
            shiny::tags$small(
                paste0("CySA ", csaVersion),
                style = paste(
                    "display:block; font-size:60%; opacity:0.75;",
                    "line-height:1.2; white-space:nowrap;"
                )
            )
        )
    )
}
