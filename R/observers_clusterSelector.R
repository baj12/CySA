# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# observers_clusterSelector.R ----
# Shiny observers for the clusterSelector app. Separating observers from
# outputs makes the server wiring easier to read and aligns with Bioconductor's
# Shiny app packaging guidelines. All observers read and write the per-session
# state environment `st` built in server_clusterSelector.R.

#' Axis Pair for a Ridge-Clicked Marker
#'
#' Decides how the main SOM 2D plot should be re-axed when the user clicks a
#' cell at \code{marker} in the ridge plot: \code{NULL} when the marker is
#' already on an axis (or unknown), the first preset pair containing the
#' marker otherwise, and as a last resort the marker on x with the first
#' other marker of \code{markers} on y. Pairs that do not repeat the marker
#' on both axes are usable; degenerate pairs are skipped.
#'
#' @param marker Clicked marker name, or \code{NULL}.
#' @param dims Current axis pair of the 2D plot.
#' @param dList Preset marker pairs.
#' @param markers All available markers.
#'
#' @return A list with \code{pair} (character length 2) and \code{idx}
#'   (position in \code{dList}, or \code{NA_integer_} for the fallback), or
#'   \code{NULL} when the axes should stay as they are.
#'
#' @keywords internal
.ridgeAxisForMarker <- function(marker, dims, dList, markers) {
    if (is.null(marker) || !marker %in% markers) {
        return(NULL)
    }
    if (marker %in% dims) {
        return(NULL)
    }
    if (length(dList) > 0L) {
        for (i in seq_along(dList)) {
            pair <- as.character(dList[[i]])
            if (length(pair) == 2L && marker %in% pair &&
                pair[1L] != pair[2L]) {
                return(list(pair = pair, idx = as.integer(i)))
            }
        }
    }
    y <- markers[!markers %in% marker][1L]
    if (is.na(y)) {
        return(NULL)
    }
    list(pair = c(marker, y), idx = NA_integer_)
}

#' Re-Axis the Main SOM 2D Plot From a Ridge-Plot Click
#'
#' Clicking an outlier dot at marker \code{m} in the ridge plot moves the
#' main 2D plot to a view that contains \code{m}: the first preset pair
#' containing it, or \code{m} on x with the first other marker on y. Clicks
#' whose payload names no marker (the "clicked cell" overlay, empty events)
#' and clicks on a marker that is already displayed leave the axes alone.
#' A matched preset pair also becomes the active preset (dropdown + static
#' views); the fallback pair behaves like a manual axis change.
#'
#' @inheritParams INTERNAL_registerClusterSelectorObservers
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.buildRidgeAxisObserver <- function(input, session, st) {
    shiny::observeEvent(
        .safeEventData("plotly_click", source = "ridgePlot",
            verbose = st$verbose),
        {
            d <- .safeEventData("plotly_click", source = "ridgePlot",
                verbose = st$verbose)
            click <- .parseRidgeClick(d)
            if (is.null(click) || is.null(click$marker)) {
                return(invisible(NULL))
            }
            dims <- c(input$currentDimX, input$currentDimY)
            move <- .ridgeAxisForMarker(click$marker, dims, st$dListRV(),
                st$rnSCE)
            if (is.null(move)) {
                return(invisible(NULL))
            }
            shiny::updateSelectInput(session = session,
                inputId = "currentDimX", selected = move$pair[1L])
            shiny::updateSelectInput(session = session,
                inputId = "currentDimY", selected = move$pair[2L])
            if (!is.na(move$idx)) {
                st$activePlot(move$idx)
                shiny::updateSelectInput(session = session,
                    inputId = "dimPairSelect",
                    selected = as.character(move$idx))
            }
        },
        ignoreNULL = TRUE
    )
    invisible(NULL)
}

#' Register clusterSelector Observers
#'
#' Attaches all `shiny::observe*` calls used by the clusterSelector app.
#' Output renderers live in `.registerClusterSelectorOutputs()`.
#'
#' @param input,output,session Shiny objects.
#' @param st Server state environment (see `.initServerState()`).
#'
#' @return Called for its side effects.
#'
#' @keywords internal
#' @name INTERNAL_registerClusterSelectorObservers
.registerClusterSelectorObservers <- function(input, output, session, st) {
    .registerStateObservers(input, session, st)
    .registerGroupObservers(input, session, st)
    .registerSampleGroupObservers(input, session, st)
    .registerPlotSelectionObservers(input, output, st)
    .buildRidgeAxisObserver(input, session, st)
    .registerDiffTestObservers(input, output, session, st)
    .registerBoxInfoObservers(input, session)
    # Close the app and return the current outputList.
    shiny::observeEvent(input$close, {
        shinyjs::js$closeWindow()
        shiny::stopApp(st$rv$outputList)
    })
    invisible(NULL)
}

#' Register the Selection-State Observers
#'
#' Keeps \code{env$outputList}, the "selected" group, the typed node ids,
#' the preset axis pair and the sample-filtered scatter data in sync with
#' the current selection.
#'
#' @inheritParams INTERNAL_registerClusterSelectorObservers
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerStateObservers <- function(input, session, st) {
    # Keep external env in sync for callers that still read env$outputList.
    shiny::observe(assign("outputList", st$rv$outputList, envir = st$env))
    # Force redraw when the implicit "selected" group is used for colouring.
    shiny::observe({
        st$rsUsed_d()
        if ("selected" %in% input$colorbyGroups) {
            shiny::isolate(st$selectedUpdate2(st$selectedUpdate2() + 1))
        }
    })
    shiny::observeEvent(input$dimPairSelect, {
        val <- input$dimPairSelect
        dList <- st$dListRV()
        if (!is.null(val) && val %in% seq_along(dList)) {
            pair <- dList[[as.integer(val)]]
            shiny::updateSelectInput(session = session,
                inputId = "currentDimX", selected = pair[1L])
            shiny::updateSelectInput(session = session,
                inputId = "currentDimY", selected = pair[2L])
            st$activePlot(as.integer(val))
        }
    })
    # Typed cluster ids update the current selection.
    shiny::observe({
        ic <- .asIntegerQuiet(st$inputClusterNumber())
        valid <- as.integer(intersect(
            as.character(ic[!is.na(ic)]), colnames(st$clusterPatientTable)
        ))
        if (length(valid) > 0L) shiny::isolate(st$rsUsed(valid))
    })
    shiny::observe(.syncSelectedGroup(st))
    # Update the subsampled data frame when the sample filter changes.
    shiny::observe({
        sampleIds <- st$sample2PlotDb()
        if (is.null(sampleIds)) {
            return(NULL)
        }
        sub <- st$sce_subsampled
        selected <- sub[, sub$sample_id %in% sampleIds]
        st$dfPlot(data.frame(t(SummarizedExperiment::assays(selected)[[1]])))
    })
    invisible(NULL)
}

#' Store the Current Selection as the "selected" Group
#'
#' Writes the current node selection into \code{outputList$selected} and
#' mirrors it in the \code{clusterNumbers} text input.
#'
#' @param st Server state environment.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.syncSelectedGroup <- function(st) {
    rs <- as.integer(unlist(st$rsUsed()))
    rs <- sort(rs[!is.na(rs)])
    outputList <- st$rv$outputList
    updateOL <- !"selected" %in% names(outputList)
    outputList$selected <- rs
    st$rv$outputList <- outputList
    if (updateOL) st$updatedoutputList()
    shiny::updateTextInput(
        inputId = "clusterNumbers",
        value = paste(rs, collapse = ", ")
    )
}

#' Register the Named-Group Observers
#'
#' Saving the selection under a name, removing one or several groups and
#' activating saved groups.
#'
#' @inheritParams INTERNAL_registerClusterSelectorObservers
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerGroupObservers <- function(input, session, st) {
    lev <- function() levels(st$sce$cluster_id)
    setGroups <- function(outputList) {
        st$rv$outputList <- .rebuildOutputList(outputList, lev())
        st$updatedoutputList()
    }
    shiny::observeEvent(input$applyName, {
        cName <- input$clusterName
        if (is.null(cName) || !nzchar(cName)) {
            return(NULL)
        }
        outputList <- st$rv$outputList
        # Typed ids are free-form text; a stale/out-of-range id saved into
        # a named group would later crash every output that indexes
        # clusterPatientTable/SOM codes by it.
        outputList[[cName]] <- .validNodeCols(st$inputClusterNumber(),
            st$clusterPatientTable)
        setGroups(outputList)
        shiny::updateSelectInput(
            inputId = "colorbyGroups",
            choices = names(st$rv$outputList),
            selected = c(input$colorbyGroups, cName)
        )
    })
    shiny::observeEvent(input$rmGrp, {
        cName <- input$clusterNameRM
        if (is.null(cName) || !nzchar(cName) ||
            !cName %in% names(st$rv$outputList)) {
            return(NULL)
        }
        outputList <- st$rv$outputList
        outputList[[cName]] <- NULL
        setGroups(outputList)
    })
    shiny::observeEvent(input$rmGroups, {
        rmGroups <- input$groupRM
        if (is.null(rmGroups) || length(rmGroups) == 0L) {
            return(NULL)
        }
        rs <- shiny::isolate(st$rsUsed())
        rmCluster <- unique(unlist(st$rv$outputList[rmGroups]))
        if (!is.null(rs)) st$rsUsed(setdiff(rs, rmCluster))
        outputList <- st$rv$outputList
        outputList[rmGroups] <- NULL
        setGroups(outputList)
    })
    .observeGroupActivation(input, session, st)
    invisible(NULL)
}

#' Activate Saved Groups in the Node-Id Input
#'
#' @inheritParams INTERNAL_registerClusterSelectorObservers
#'
#' @return The observer, invisibly.
#'
#' @keywords internal
.observeGroupActivation <- function(input, session, st) {
    shiny::observeEvent(input$clusterNameSelect, {
        outputList <- st$rv$outputList
        combined <- unique(unlist(
            outputList[names(outputList) %in% input$clusterNameSelect]
        ))
        if (length(combined) > 0L) {
            shiny::updateTextInput(session = session,
                inputId = "clusterNumbers",
                value = paste(combined, collapse = ", "))
        }
    })
}

#' Register the t-Test Group Observers
#'
#' Updates the \code{group1} / \code{group2} choices when the grouping
#' variable changes and keeps the two groups mutually exclusive.
#'
#' @inheritParams INTERNAL_registerClusterSelectorObservers
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerSampleGroupObservers <- function(input, session, st) {
    levs <- function() {
        groupsVar <- input$groupsVar
        info <- st$metaD$experiment_info
        if (is.null(groupsVar) || !groupsVar %in% colnames(info)) {
            return(NULL)
        }
        levels(info[, groupsVar])
    }
    shiny::observeEvent(input$groupsVar, {
        l <- levs()
        if (is.null(l)) {
            return(NULL)
        }
        shiny::updateSelectInput(session = session, inputId = "group1",
            choices = l)
        shiny::updateSelectInput(session = session, inputId = "group2",
            choices = l)
    })
    # Each group's choices exclude the levels picked in the other group.
    exclude <- function(from, to) {
        shiny::observe({
            l <- levs()
            if (is.null(l)) {
                return(NULL)
            }
            picked <- input[[from]]
            keep <- shiny::isolate(input[[to]])
            shiny::updateSelectInput(session = session, inputId = to,
                choices = l[!l %in% picked], selected = keep)
        })
    }
    exclude("group1", "group2")
    exclude("group2", "group1")
    invisible(NULL)
}

#' Register the Plot Selection and Zoom Observers
#'
#' Lasso / box selections and clicks on the SOM, dimension-reduction,
#' raster, dendrogram, scatter and FlowSOM star plots update the node
#' selection according to \code{input$selectMode}; relayout events of the
#' SOM plots update the zoom state.
#'
#' @inheritParams INTERNAL_registerClusterSelectorObservers
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerPlotSelectionObservers <- function(input, output, st) {
    verbose <- st$verbose
    .buildSOMDataObservers(st$nPlots, input, output, st$rsUsed,
        st$rsUsed_d, st$activePlot, .inputSelect, verbose)
    .buildZoomObservers(st$nPlots, input, st$zoomFunc, verbose)
    shiny::observe(st$zoomFunc(.safeEventData("plotly_relayout",
        source = "somDataMain", verbose = verbose), 1L))
    for (src in c("somDataMain", "tsne", "umap", "pca", "nodeProps")) {
        .buildSelectionObserver(src, input, st$rsUsed, st$rsUsed_d,
            .inputSelect, verbose)
    }
    # The raster grid selects relative to the undebounced selection.
    .buildSelectionObserver("somGrid", input, st$rsUsed, st$rsUsed,
        .inputSelect, verbose)
    sources <- c("dendPlotly", if (!is.null(st$fsom)) "flowSOMStars")
    for (src in sources) {
        for (event in c("plotly_selected", "plotly_click")) {
            .observeCustomdataSelection(input, st, src, event)
        }
    }
    shiny::observeEvent(
        .safeEventData("plotly_selected", source = "scatterPlot",
            verbose = verbose),
        {
            d <- .safeEventData("plotly_selected", source = "scatterPlot",
                verbose = verbose)
            newRs <- .scatterSelectionToNodes(st, input, d)
            if (!is.null(newRs)) shiny::isolate(st$rsUsed(newRs))
        }
    )
    invisible(NULL)
}

#' Observe a Plotly Event Carrying Node Ids in \code{customdata}
#'
#' @param input Shiny input object.
#' @param st Server state environment.
#' @param source Plotly event source id.
#' @param event Plotly event name, e.g. \code{"plotly_click"}.
#'
#' @return The observer, invisibly.
#'
#' @keywords internal
.observeCustomdataSelection <- function(input, st, source, event) {
    force(source)
    force(event)
    shiny::observeEvent(
        .safeEventData(event, source = source, verbose = st$verbose),
        {
            d <- .safeEventData(event, source = source,
                verbose = st$verbose)
            if (is.null(d) || NROW(d) == 0L || is.null(d$customdata)) {
                return(NULL)
            }
            newRs <- .applySelectMode(
                shiny::isolate(st$rsUsed_d()),
                as.integer(d$customdata),
                shiny::isolate(input$selectMode)
            )
            shiny::isolate(st$rsUsed(newRs))
        }
    )
}

#' SOM Nodes of the Cells Inside a Scatter-Plot Box Selection
#'
#' @param st Server state environment.
#' @param input Shiny input object.
#' @param d \code{plotly_selected} event data of the scatter plot.
#'
#' @return The new node selection after applying \code{input$selectMode},
#'   or \code{NULL} when the selection is empty or has zero area.
#'
#' @keywords internal
.scatterSelectionToNodes <- function(st, input, d) {
    rs <- shiny::isolate(st$rsUsed_d())
    dimSel <- st$dimSelection()
    sampleIds <- shiny::isolate(input$samples2plot)
    plotIdx <- st$activePlot()
    if (is.null(d)) {
        return(NULL)
    }
    shiny::req(d$curveNumber)
    d <- d[d$curveNumber == 1L, , drop = FALSE]
    # Zero-area (degenerate) selections from a single point or colinear
    # points are ignored rather than wiping the current selection.
    if (nrow(d) == 0L || min(d$x) == max(d$x) || min(d$y) == max(d$y)) {
        return(NULL)
    }
    dfLocal <- st$dfPlot()
    dims <- make.names(dimSel[[plotIdx]]$dims)
    xs <- dfLocal[, dims[1L], drop = TRUE]
    ys <- dfLocal[, dims[2L], drop = TRUE]
    ids <- which(xs > min(d$x) & xs < max(d$x) &
        ys > min(d$y) & ys < max(d$y))
    sub <- st$sce_subsampled
    ids <- SingleCellExperiment::colData(
        sub[, sub$sample_id %in% sampleIds]
    )[ids, "cluster_id"]
    .applySelectMode(rs, unique(as.integer(ids)),
        shiny::isolate(input$selectMode))
}
