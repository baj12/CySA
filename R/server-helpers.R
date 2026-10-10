# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server-helpers.R ----
# Non-reactive helper functions used inside the clusterSelector server. Moving
# these out of the server closure reduces the size of clusterSelector() and
# makes each piece independently testable.

#' Synchronize Select Input Choices with outputList
#'
#' Updates the UI select inputs that depend on the names of \code{outputList}.
#'
#' @param session Shiny session object.
#' @param input Shiny input object.
#' @param outputList Named list of cluster groupings.
#' @param metaD \code{metadata(sce)} list.
#'
#' @keywords internal
.updateOutputListInputs <- function(session, input, outputList, metaD) {
    olNames <- names(outputList)
    update <- function(id, choices, selected = isolate(input[[id]])) {
        shiny::updateSelectInput(session = session, id, choices = choices,
            selected = selected)
    }
    update("clusterNameSelect", olNames)
    numCols <- unlist(lapply(metaD$experiment_info, is.numeric),
        use.names = FALSE)
    eI <- apply(metaD$experiment_info[, numCols, drop = FALSE], 2, as.numeric)
    if (any(is.na(eI))) {
        stop("NAs produced when converting experiment_info to numeric")
    }
    choices <- c("none", colnames(eI), olNames)
    update("compareStatsTo", choices)
    update("relativeTo", choices)
    newSel <- union(isolate(input$upsetSelection),
        intersect(c("Rest", "selected"), olNames))
    if (length(newSel) == 0) newSel <- olNames
    update("upsetSelection", olNames, newSel)
    update("colorbyGroups", olNames)
    update("groupRM", olNames)
    shiny::updateSelectInput(session = session, "clusterNameRM",
        choices = olNames)
    update("diagGroups", olNames)
    update("orderedHeatmapGroups", olNames)
    update("diffGroups", olNames)
}


#' Compute Relative Cell Counts
#'
#' Computes the percentage of selected cells relative to a reference column,
#' another cluster group, or the total population.
#'
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param rs Selected cluster ids.
#' @param relativeToCol Reference column or group name.
#' @param expInfo Sample-level metadata data.frame.
#' @param outputList Named list of cluster groupings.
#'
#' @return Numeric vector of percentages.
#'
#' @keywords internal
.computeRelativeCounts <- function(clusterPatientTable, rs, relativeToCol,
    expInfo, outputList) {
    rs <- .validNodeCols(rs, clusterPatientTable)
    switch(relativeToCol,
        "none" = rowSums(clusterPatientTable[, rs, drop = FALSE]) /
            rowSums(clusterPatientTable[, , drop = FALSE]) * 100,
        {
            if (relativeToCol %in% colnames(expInfo)) {
                rowSums(clusterPatientTable[, rs, drop = FALSE]) /
                    expInfo[rownames(clusterPatientTable), relativeToCol] * 100
            } else if (relativeToCol %in% names(outputList)) {
                grpNodes <- .validNodeCols(outputList[[relativeToCol]],
                    clusterPatientTable)
                rowSums(clusterPatientTable[, rs, drop = FALSE]) /
                    rowSums(clusterPatientTable[, grpNodes,
                        drop = FALSE]) * 100
            } else {
                stop("relativeToCol '", relativeToCol, "' not found")
            }
        }
    )
}

#' Node Ids That Are Valid Columns of clusterPatientTable
#'
#' Saved groups (\code{outputList} entries) can go stale relative to the
#' current data -- e.g. a node id typed and named before the selection was
#' revalidated, or an \code{outputList} supplied to \code{clusterSelector()}
#' from an earlier session with a different SOM. Indexing
#' \code{clusterPatientTable} by such an id throws "subscript out of
#' bounds" and takes the whole output down with it, so every caller that
#' indexes by external node ids filters them through this first.
#'
#' @param nodes Requested node ids.
#' @param clusterPatientTable Table of sample by cluster counts.
#'
#' @return \code{nodes}, restricted to ids that are actual columns of
#'   \code{clusterPatientTable}.
#'
#' @keywords internal
.validNodeCols <- function(nodes, clusterPatientTable) {
    as.integer(intersect(as.character(nodes),
        colnames(clusterPatientTable)))
}


#' Update outputList After Cluster Assignment Change
#'
#' Rebuilds \code{outputList} after a group is removed or renamed, ensuring the
#' \code{"Rest"} group is always correct.
#'
#' @param outputList Current named list of cluster groupings.
#' @param sceLevels Character or integer vector of all cluster levels.
#' @param removed Optional vector of cluster ids that should be excluded.
#'
#' @return Updated named list.
#'
#' @keywords internal
.rebuildOutputList <- function(outputList, sceLevels, removed = NULL) {
    if (!is.null(removed)) {
        for (nm in names(outputList)) {
            outputList[[nm]] <- setdiff(outputList[[nm]], removed)
        }
    }

    # The previous "Rest" must not count as used, or Rest would always end
    # up empty (and be dropped) after the first group change.
    used <- unique(unlist(outputList[names(outputList) != "Rest"]))
    outputList[["Rest"]] <- as.integer(sceLevels[!sceLevels %in% used])

    for (na in names(outputList)) {
        if (length(outputList[[na]]) == 0) outputList[[na]] <- NULL
    }
    outputList
}

#' Combine an existing selection with a newly picked set of IDs
#'
#' Shared "selectMode" semantics used by every plotly click/select observer
#' in this app. \code{"remove others"} narrows the current selection down to
#' only the overlap with the new pick (it does NOT replace the selection
#' with the new pick alone).
#'
#' @param rs Integer vector of the currently selected IDs.
#' @param picked Integer vector of newly picked IDs.
#' @param mode One of \code{"remove others"}, \code{"add"}, \code{"remove"},
#'   or anything else (falls through to \code{picked}). \code{NULL} is
#'   treated as \code{"view"} (falls through to \code{picked}) rather than
#'   erroring in \code{switch()}.
#'
#' @return Integer vector: the combined selection.
#' @keywords internal
.applySelectMode <- function(rs, picked, mode) {
    mode <- if (is.null(mode)) "view" else mode
    # Defensive: empty picks should never modify the current selection.
    if (length(picked) == 0L) {
        return(rs)
    }
    switch(
        EXPR = mode,
        "remove others" = intersect(picked, rs),
        "add" = unique(c(rs, picked)),
        "remove" = setdiff(rs, picked),
        picked
    )
}

# inputSelect = function ----
.inputSelect <- function(d, rs, mode) {
    if (is.null(rs)) {
        return(integer(0))
    }
    if (is.null(d)) {
        return(as.integer(rs))
    }

    node_ids <- NULL

    if (is.data.frame(d)) {
        if ("key" %in% names(d) && length(d$key) > 0L && !all(is.na(d$key))) {
            node_ids <- as.integer(unlist(d$key))
        } else if ("customdata" %in% names(d) && !all(is.na(d$customdata))) {
            node_ids <- as.integer(unlist(d$customdata))
        } else if ("pointNumber" %in% names(d)) {
            if ("curveNumber" %in% names(d)) {
                d <- d[d$curveNumber == 0L, , drop = FALSE]
                if (nrow(d) == 0L) {
                    return(as.integer(rs))
                }
            }
            node_ids <- as.integer(unlist(d$pointNumber)) + 1L
        }
    } else {
        node_ids <- .asIntegerQuiet(unlist(d))
    }

    if (is.null(node_ids)) {
        return(as.integer(rs))
    }
    node_ids <- node_ids[!is.na(node_ids) & node_ids > 0L]
    if (length(node_ids) == 0L) {
        return(as.integer(rs))
    }

    as.integer(.applySelectMode(rs, node_ids, mode))
}


#' Write All ClusterSelector Plots to a PDF
#'
#' Builds the full "Download Plots" PDF export: the dendrogram, count/percent
#' bar charts, every preset SOM marker-pair panel (not just the currently
#' active one), the SOM raster overlay, t-SNE/UMAP/PCA, the main scatter plot,
#' FlowSOM pie/star plots when available, both violin plots, and the UpSet
#' plot. Each entry in \code{plotRegistry} is wrapped in its own
#' \code{tryCatch()} so one broken/unavailable plot is skipped with a warning
#' rather than aborting the whole export (which previously surfaced to users
#' as a corrupt/HTML "download").
#'
#' @param file Output PDF path (as required by \code{downloadHandler}'s
#'   \code{content} argument).
#' @param rs Integer vector of currently selected SOM node ids.
#' @param sce Full \code{SingleCellExperiment}.
#' @param sceRN,sceCN \code{rownames(sce)}, \code{colnames(sce)}.
#' @param metaD \code{metadata(sce)} list.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param dendPlotObj A \code{dendrogram} object (already computed; this
#'   function only plots it).
#' @param dimPairs List of marker-pair vectors, one per SOM panel shown on
#'   screen (i.e. \code{dListRV()}) -- ALL pairs, not just the active plot.
#' @param getBasePlotFn Function \code{(x, y, colorVar, sizeVar) -> ggplot}
#'   used to build each SOM panel's base plot, matching whatever the live
#'   static-grid renderer uses.
#' @param dlColorVar,dlSizeVar SOM statistic used for point colour / size.
#' @param pctl Central fraction of the SOM codes kept by the axis limits.
#' @param somRasterXy Data frame of selected-node raster coordinates (or
#'   \code{NULL}).
#' @param baseRasterGgplot Pre-built raster background ggplot (or
#'   \code{NULL}).
#' @param plotRegistry Named list of zero-argument functions, each returning
#'   a plot object (\code{ggplot}, base plot, or a \code{plotly}/htmlwidget
#'   whose \code{print()} method renders sensibly) or \code{NULL} to skip.
#'   Errors thrown by any entry are caught and reported as a warning, not
#'   propagated.
#'
#' @keywords internal
.writeClusterSelectorPdf <- function(file, rs, sce, sceRN, sceCN, metaD,
    somCodesName, dendPlotObj, dimPairs, getBasePlotFn, dlColorVar,
    dlSizeVar, pctl, somRasterXy, baseRasterGgplot, plotRegistry = list()) {
    grDevices::pdf(file = file)
    on.exit(grDevices::dev.off(), add = TRUE)
    .pdfPage("Skipping 'dendrogram': ",
        function() plot(dendPlotObj, main = "dendrogram"))
    somCodes <- metaD[[somCodesName]]
    for (i in seq_along(dimPairs)) {
        pair <- dimPairs[[i]]
        .pdfPage(
            paste0("Skipping SOM panel ", i, " (",
                paste(pair, collapse = "-"), "): "),
            function() {
                pp1 <- getBasePlotFn(pair[1L], pair[2L], dlColorVar, dlSizeVar)
                methods::show(ggsomPlot(pp1, 1L, rs,
                    dimSelection = list(list(dims = pair)),
                    sce = sce, metaD = metaD,
                    xlim = .quantileLimits(somCodes, pair[1L], pctl),
                    ylim = .quantileLimits(somCodes, pair[2L], pctl)))
            }
        )
    }
    if (!is.null(somRasterXy) && !is.null(baseRasterGgplot)) {
        .pdfPage("Skipping SOM raster overlay in PDF export: ", function() {
            methods::show(.rasterWithSelection(baseRasterGgplot, somRasterXy))
        })
    }
    # Registry-driven plots (tSNE/UMAP/PCA/scatter/FlowSOM/violins/UpSet).
    for (nm in names(plotRegistry)) {
        fn <- plotRegistry[[nm]]
        if (is.null(fn)) next
        p <- tryCatch(fn(), error = function(e) {
            warning("Skipping '", nm, "' in PDF export: ", conditionMessage(e))
            NULL
        })
        if (!is.null(p)) {
            .pdfPage(
                paste0("Skipping '", nm, "' in PDF export (print failed): "),
                function() methods::show(p)
            )
        }
    }
    invisible(NULL)
}

#' Draw One PDF Page, Turning Errors Into Warnings
#'
#' @param prefix Warning text placed before the error message.
#' @param draw Function drawing the page.
#'
#' @return Called for its side effect.
#'
#' @keywords internal
.pdfPage <- function(prefix, draw) {
    tryCatch(draw(), error = function(e) {
        warning(prefix, conditionMessage(e))
    })
    invisible(NULL)
}


#' Convert to Integer Without Coercion Warnings
#'
#' Like \code{as.integer()}, but values that do not look like numbers become
#' \code{NA} silently instead of raising "NAs introduced by coercion".
#'
#' @param x Vector to convert.
#'
#' @return An integer vector of the same length as \code{x}.
#'
#' @keywords internal
.asIntegerQuiet <- function(x) {
    if (is.numeric(x)) {
        return(as.integer(x))
    }
    x <- trimws(as.character(unlist(x)))
    out <- rep(NA_integer_, length(x))
    ok <- !is.na(x) &
        grepl("^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$", x)
    out[ok] <- as.integer(as.numeric(x[ok]))
    out
}
