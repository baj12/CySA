# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# diffcyt-helpers.R ----
# diffcyt differential-abundance (DA) and differential-state (DS) testing on
# CySA's own saved cluster groups (rv$outputList), rather than on raw SOM
# nodes -- CySA's equivalent of CATALYST's manually merged metaclusters.
# vignettes/catalyst_workflow.Rmd reproduces the reference diffcyt workflow
# (on CATALYST metaclusters) that this is adapted from.
#
# diffcyt's calcCounts()/calcMedians() expect a "reversed" SummarizedExperiment
# (rows = cells, columns = markers), with per-cell `sample_id`/`cluster_id` in
# rowData() and per-marker `marker_class`/`marker_name` in colData() -- the
# same convention HDCytoData/prepareData() use, not the usual SCE convention.

#' Map SOM Node Ids to Their Saved Group Name
#'
#' @param outputList Named list of SOM node ids per saved group.
#'
#' @return A character vector named by node id, valued by group name.
#'
#' @keywords internal
.nodeGroupMap <- function(outputList) {
    nodes <- unlist(lapply(outputList, as.character), use.names = FALSE)
    groups <- rep(names(outputList), lengths(outputList))
    stats::setNames(groups, nodes)
}

#' Check That the Compared Groups Do Not Share Any SOM Node
#'
#' A differential test needs each cluster/node assigned to exactly one of
#' the compared groups. CySA's saved groups can overlap (most commonly:
#' \code{"selected"} mirrors whatever nodes were just named into a new
#' group), and \code{.nodeGroupMap()}'s node-to-group lookup silently keeps
#' only the first group for a shared node -- leaving the other group with
#' zero cells, an all-\code{NA} row in
#' \code{diffcyt::calcMediansByClusterMarker()}, and a crash in
#' \code{diffcyt::plotHeatmap()}'s internal \code{hclust()} row clustering.
#' Checking up front turns that crash into a clear error.
#'
#' @param outputList Named list of SOM node ids per saved group (already
#'   restricted to the groups being compared).
#'
#' @return Invisible \code{NULL}; errors if any two groups overlap.
#'
#' @keywords internal
.checkGroupsDisjoint <- function(outputList) {
    pairs <- utils::combn(names(outputList), 2, simplify = FALSE)
    problems <- Filter(Negate(is.null), lapply(pairs, function(p) {
        shared <- intersect(as.character(outputList[[p[1]]]),
            as.character(outputList[[p[2]]]))
        if (length(shared) == 0L) {
            return(NULL)
        }
        paste0(p[1], " & ", p[2], " share node(s) ",
            paste(shared, collapse = ", "))
    }))
    if (length(problems) > 0L) {
        stop(
            "The compared groups overlap, so some SOM nodes do not belong ",
            "to exactly one group (each node must belong to only one of ",
            "the groups being compared): ",
            paste(problems, collapse = "; "), ". A common cause is ",
            "leaving 'selected' selected alongside the group you just ",
            "named from it -- deselect it or any other overlapping group."
        )
    }
    invisible(NULL)
}

#' Saved Groups That Actually Contain at Least One Cell
#'
#' A saved group with zero matching cells (e.g. an emptied-out "Rest", or a
#' group created but never populated) produces an all-\code{NA} row in
#' \code{diffcyt::calcMediansByClusterMarker()}'s output, which crashes
#' \code{diffcyt::plotHeatmap()}'s internal row clustering
#' (\code{hclust(): "NA/NaN/Inf in foreign function call"}). Dropping empty
#' groups before building the diffcyt object avoids this entirely.
#'
#' @param sce Full \code{SingleCellExperiment} (colData \code{cluster_id}).
#' @param outputList Named list of SOM node ids per saved group.
#'
#' @return \code{outputList}, restricted to groups with at least one cell.
#'
#' @keywords internal
.nonEmptyGroups <- function(sce, outputList) {
    nodeIds <- as.character(SingleCellExperiment::colData(sce)$cluster_id)
    present <- vapply(outputList, function(nodes) {
        any(as.character(nodes) %in% nodeIds)
    }, logical(1))
    outputList[present]
}

#' Build a diffcyt-Shaped SummarizedExperiment From the Saved Groups
#'
#' \code{diffcyt::plotHeatmap()} always draws a cluster-identity panel from
#' the \code{"type"}-class markers, even for a DA test, so \code{typeMarkers}
#' (normally the SOM clustering markers) are always carried along;
#' \code{stateMarkers} (the markers a DS test is run on) are additional and
#' only relevant for DS.
#'
#' @param sce Full \code{SingleCellExperiment} (colData \code{sample_id},
#'   \code{cluster_id}).
#' @param outputList Named list of SOM node ids per saved group; each saved
#'   group becomes one diffcyt "cluster".
#' @param experimentInfo \code{metadata(sce)$experiment_info}, with a
#'   \code{sample_id} column.
#' @param typeMarkers Marker (row) names carried as \code{marker_class ==
#'   "type"} (cluster-identity markers).
#' @param stateMarkers Marker (row) names carried as \code{marker_class ==
#'   "state"} (markers tested in a DS test); \code{character(0)} for DA.
#'
#' @return A \code{SummarizedExperiment} in diffcyt's own layout.
#'
#' @keywords internal
.buildDiffcytSce <- function(sce, outputList, experimentInfo, typeMarkers,
    stateMarkers = character(0)) {
    nodeToGroup <- .nodeGroupMap(outputList)
    cd <- SingleCellExperiment::colData(sce)
    groupId <- nodeToGroup[as.character(cd$cluster_id)]
    keep <- !is.na(groupId)
    if (!any(keep)) {
        stop("None of the cells belong to a saved group; save at least ",
            "one group before running a differential test.")
    }
    markers <- union(typeMarkers, stateMarkers)
    exprMat <- t(SummarizedExperiment::assay(sce, 1)[markers, keep,
        drop = FALSE])
    rowD <- S4Vectors::DataFrame(
        sample_id = factor(cd$sample_id[keep],
            levels = as.character(experimentInfo$sample_id)),
        cluster_id = factor(groupId[keep], levels = names(outputList))
    )
    colD <- S4Vectors::DataFrame(
        marker_name = markers,
        marker_class = factor(
            ifelse(markers %in% stateMarkers, "state", "type"),
            levels = c("type", "state", "none")
        ),
        row.names = markers
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = exprMat), rowData = rowD, colData = colD,
        metadata = list(experiment_info = experimentInfo)
    )
}

#' Run a diffcyt Differential Test on the Saved Groups
#'
#' Differential abundance (DA, GLMM) or differential state (DS, LMM) test of
#' \code{outputList}'s saved groups across samples, following diffcyt's own
#' low-level pipeline (\code{calcCounts()}/\code{calcMedians()} +
#' \code{testDA_GLMM()}/\code{testDS_LMM()}) rather than its CATALYST-specific
#' wrapper, since the clustering here is already fixed by the user's own
#' curated groups.
#'
#' @param sce Full \code{SingleCellExperiment}.
#' @param outputList Named list of SOM node ids per saved group.
#' @param experimentInfo \code{metadata(sce)$experiment_info}.
#' @param fixedCol Two-level factor column of \code{experimentInfo} used as
#'   the fixed effect (the comparison of interest).
#' @param randomCol Optional factor column of \code{experimentInfo} added as
#'   a random intercept, or \code{NULL}/\code{"none"} for no random effect.
#' @param analysisType \code{"DA"} or \code{"DS"}.
#' @param markers Markers to test (\code{DS} only; ignored for \code{DA},
#'   which tests cluster abundance, not marker expression). Defaults to all
#'   of \code{rownames(sce)}.
#' @param typeMarkers Cluster-identity markers carried along for
#'   \code{diffcyt::plotHeatmap()}'s cell-type panel (used for both DA and
#'   DS). Defaults to \code{metadata(sce)$map$colsUsed}, falling back to
#'   \code{rownames(sce)}.
#'
#' @return A list with \code{res}, \code{d_counts}, \code{d_medians} (DS
#'   only), \code{d_medians_by_cluster_marker} and \code{analysis_type}, fed
#'   to \code{\link{diffTestTable}()}/\code{\link{plotDiffTestHeatmap}()}.
#'
#' @keywords internal
runDiffTest <- function(sce, outputList, experimentInfo, fixedCol,
    randomCol = NULL, analysisType = c("DA", "DS"), markers = NULL,
    typeMarkers = NULL) {
    analysisType <- match.arg(analysisType)
    nonEmpty <- .nonEmptyGroups(sce, outputList)
    if (length(nonEmpty) < length(outputList)) {
        warning(
            "Dropped empty saved group(s) with no matching cells: ",
            paste(setdiff(names(outputList), names(nonEmpty)),
                collapse = ", "),
            call. = FALSE
        )
    }
    if (length(nonEmpty) < 2L) {
        stop("At least 2 non-empty saved groups are required for a ",
            "differential test; found ", length(nonEmpty), ".")
    }
    outputList <- nonEmpty
    .checkGroupsDisjoint(outputList)
    if (is.null(typeMarkers)) {
        typeMarkers <- S4Vectors::metadata(sce)$map$colsUsed %||%
            rownames(sce)
    }
    stateMarkers <- if (analysisType == "DS") {
        if (is.null(markers)) rownames(sce) else markers
    } else {
        character(0)
    }
    # diffcyt's calcCounts()/calcMedians()/topTable() sort and subset
    # samples via levels(colData(d_counts)$sample_id); experiment_info's
    # sample_id must therefore be a factor, not the character column
    # CySA_example_sce()/HDCytoData-style metadata typically carries.
    experimentInfo$sample_id <- factor(as.character(experimentInfo$sample_id))
    colsRandom <- if (is.null(randomCol) || identical(randomCol, "none")) {
        NULL
    } else {
        randomCol
    }
    dSe <- .buildDiffcytSce(sce, outputList, experimentInfo, typeMarkers,
        stateMarkers)
    dCounts <- diffcyt::calcCounts(dSe)
    frm <- diffcyt::createFormula(experimentInfo, cols_fixed = fixedCol,
        cols_random = colsRandom)
    contrast <- diffcyt::createContrast(c(0, 1))
    if (analysisType == "DA") {
        res <- diffcyt::testDA_GLMM(dCounts, frm, contrast)
        dMedians <- NULL
    } else {
        dMedians <- diffcyt::calcMedians(dSe)
        res <- diffcyt::testDS_LMM(dCounts, dMedians, frm, contrast)
    }
    list(
        res = res, d_counts = dCounts, d_medians = dMedians,
        d_medians_by_cluster_marker =
            diffcyt::calcMediansByClusterMarker(dSe),
        analysis_type = analysisType
    )
}

#' Top-Table of a Differential Test Result
#'
#' @param diffResult Result of \code{\link{runDiffTest}()}.
#' @param top_n Number of top (most significant) rows to keep.
#'
#' @return A data frame, as returned by \code{diffcyt::topTable()}.
#'
#' @keywords internal
diffTestTable <- function(diffResult, top_n = 20) {
    diffcyt::topTable(
        res = diffResult$res, d_counts = diffResult$d_counts,
        d_medians = diffResult$d_medians,
        order_by = "p_adj", top_n = top_n,
        show_props = identical(diffResult$analysis_type, "DA"),
        show_meds = identical(diffResult$analysis_type, "DS"),
        format_vals = TRUE, digits = 3
    )
}

#' Heatmap of a Differential Test Result
#'
#' @inheritParams diffTestTable
#' @param threshold Adjusted p-value threshold used to flag rows.
#'
#' @return A \code{ComplexHeatmap} object (drawn via \code{grid} graphics; can
#'   be passed directly to \code{shiny::renderPlot()}, as CySA's existing
#'   UpSet plot already does).
#'
#' @keywords internal
plotDiffTestHeatmap <- function(diffResult, top_n = 20, threshold = 0.1) {
    diffcyt::plotHeatmap(
        res = diffResult$res, analysis_type = diffResult$analysis_type,
        d_counts = diffResult$d_counts, d_medians = diffResult$d_medians,
        d_medians_by_cluster_marker = diffResult$d_medians_by_cluster_marker,
        top_n = top_n, threshold = threshold
    )
}
