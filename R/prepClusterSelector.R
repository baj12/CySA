# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

#' Default Marker Pairs for 2D Plots
#'
#' Builds a default list of marker pairs from a \code{SingleCellExperiment}
#' object.
#'
#' @param sce A \code{SingleCellExperiment}.
#'
#' @return A list of character pairs.
#'
#' @keywords internal
.defaultDList <- function(sce) {
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed
    rn <- if (!is.null(colsUsed) && length(colsUsed) >= 12) {
        colsUsed
    } else {
        rownames(sce)
    }
    if (length(rn) < 12) {
        stop("sce must have at least 12 row names to build default dList")
    }
    list(
        d1 = c(rn[1], rn[2]),
        d2 = c(rn[3], rn[4]),
        d3 = c(rn[5], rn[6]),
        d4 = c(rn[7], rn[8]),
        d5 = c(rn[9], rn[10]),
        d6 = c(rn[11], rn[12])
    )
}


#' Prepare Data for the Cluster Selector Shiny App
#'
#' Subsamples a \code{SingleCellExperiment} object and builds the inputs
#' required by \code{\link{clusterSelector}}.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing the full dataset.
#' @param somFile Path to the SOM object. Used to cache/restore the subsampled
#'   version.
#' @param dList Optional list of marker pairs for 2D plots. If \code{NULL},
#'   defaults to the first 12 row names of \code{sce}.
#' @param total_cells_to_sample Number of cells to subsample in total.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with \code{sce}, \code{sce_subsampled}, and \code{dList}.
#'
#' @examples
#' sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
#' prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
#' names(prepped)
#'
#' @export
prepClusterSelectorData <- function(sce, somFile = NULL, dList = NULL,
    total_cells_to_sample = 100000, somCodesName = "SOM_codes",
    assay = "exprs", seed = 123) {
    if (length(total_cells_to_sample) != 1 ||
        !is.numeric(total_cells_to_sample) || total_cells_to_sample < 1) {
        stop("'total_cells_to_sample' must be a positive scalar")
    }
    if (is.null(dList)) {
        dList <- .defaultDList(sce)
    }
    cacheFile <- if (!is.null(somFile)) {
        paste0(tools::file_path_sans_ext(somFile), ".subsampled.RData")
    }
    if (!is.null(cacheFile) && file.exists(cacheFile)) {
        env <- new.env()
        load(cacheFile, envir = env)
        sce_subsampled <- env$sce_subsampled
    } else {
        sce_subsampled <- withr::with_seed(seed, .subsampleSce(
            sce, total_cells_to_sample, somCodesName
        ))
        if (!is.null(cacheFile)) {
            save(sce_subsampled, file = cacheFile)
        }
    }
    sce_subsampled <- .restoreSubsampledAssay(sce_subsampled, sce, assay)
    list(sce = sce, sce_subsampled = sce_subsampled, dList = dList)
}

#' Stratified Subsample of the Cells
#'
#' Draws from every (sample, cluster) group proportionally to its size, plus
#' ten cells, so that small groups stay represented.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param total Approximate total number of cells to draw.
#' @param somCodesName Name of the SOM codes metadata slot.
#'
#' @return The subsampled \code{SingleCellExperiment}; \code{cluster_id} is
#'   a factor with one level per SOM node.
#'
#' @keywords internal
.subsampleSce <- function(sce, total, somCodesName) {
    cd <- as.data.frame(SingleCellExperiment::colData(sce))
    groups <- dplyr::summarise(
        dplyr::group_by(cd, .data$sample_id, .data$cluster_id),
        group_size = dplyr::n(), .groups = "drop"
    )
    nToSample <- ceiling(groups$group_size / sum(groups$group_size) * total +
        10)
    parts <- lapply(seq_len(nrow(groups)), function(g) {
        size <- groups$group_size[g]
        idx <- if (nToSample[g] > size) seq_len(size) else sample(size,
            nToSample[g])
        inGroup <- cd$sample_id == as.character(groups$sample_id[g]) &
            cd$cluster_id == as.character(groups$cluster_id[g])
        sce[, inGroup][, idx]
    })
    sub <- do.call(SummarizedExperiment::cbind, parts)
    # Downstream helpers need a factor cluster_id with all node levels.
    sub$cluster_id <- factor(
        as.character(SingleCellExperiment::colData(sub)$cluster_id),
        levels = seq_len(nrow(S4Vectors::metadata(sce)[[somCodesName]]))
    )
    sub
}

#' Restore the Assay of a Subsampled SCE Saved Without Assays
#'
#' @param sub Subsampled \code{SingleCellExperiment} (with \code{id} in
#'   \code{colData} indexing the full object).
#' @param sce The full \code{SingleCellExperiment}.
#' @param assay Name for the restored assay.
#'
#' @return \code{sub}, with the first assay of \code{sce} restored when it
#'   had none.
#'
#' @keywords internal
.restoreSubsampledAssay <- function(sub, sce, assay) {
    if (!is.null(sub) && length(SummarizedExperiment::assays(sub)) == 0) {
        SummarizedExperiment::assays(sub)[[1]] <-
            SummarizedExperiment::assays(sce)[[1]][,
                SummarizedExperiment::colData(sub)$id]
        names(SummarizedExperiment::assays(sub)) <- assay
    }
    sub
}


#' Build SOM Codes From Cell-Level Expression Data
#'
#' Computes per-cluster mean expression for all markers in \code{sce}.
#' Markers already present in the existing SOM codes metadata slot are left
#' unchanged; missing markers are calculated from the cell-level assay.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id} in
#'   \code{colData} and SOM codes stored in
#'   \code{metadata(sce)[[somCodesName]]}.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param markerList Optional character vector of markers to include. If
#'   \code{NULL}, all markers except \code{"label"} and \code{"TIME"} are used.
#' @param assay Name of the assay to use.
#'
#' @return The updated \code{SingleCellExperiment} with refreshed SOM codes.
#'
#' @keywords internal
.buildSOMCodes <- function(sce,
    somCodesName = "SOM_codes",
    markerList = NULL,
    assay = "exprs") {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]

    if (is.null(markerList)) {
        markerList <- setdiff(rownames(sce), c("label", "TIME"))
    }

    # Only recompute markers that are not already present in SOM codes.
    markers_to_compute <- setdiff(markerList, colnames(somCodes))
    if (length(markers_to_compute) == 0) {
        return(sce)
    }

    exprs <- SummarizedExperiment::assay(sce, assay)[markers_to_compute, ,
        drop = FALSE]
    cluster_id <- as.integer(SingleCellExperiment::colData(sce)$cluster_id)
    n_clusters <- nrow(somCodes)

    # Sum expression per cluster, then divide by cluster counts.
    cluster_sums <- rowsum(t(exprs), cluster_id, na.rm = TRUE)
    cluster_counts <- tabulate(cluster_id, nbins = n_clusters)
    # Avoid division by zero; the sums are 0 anyway.
    cluster_counts[cluster_counts == 0] <- 1
    cluster_means <- cluster_sums / cluster_counts

    # rowsum may drop empty clusters; ensure full n_clusters rows.
    if (nrow(cluster_means) < n_clusters) {
        full <- matrix(0, nrow = n_clusters, ncol = length(markers_to_compute))
        colnames(full) <- markers_to_compute
        present <- as.integer(rownames(cluster_means))
        full[present, ] <- as.matrix(cluster_means)
        cluster_means <- full
    }

    new_codes <- cbind(somCodes, cluster_means)
    # Reorder to the requested marker order.
    S4Vectors::metadata(sce)[[somCodesName]] <-
        new_codes[, markerList, drop = FALSE]
    sce
}


#' Build SOM Summary Statistics
#'
#' Computes per-cluster summary statistics (median, mean, third quartile, max,
#' cell count) based on distances between SOM codes and their assigned cells.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id} in
#'   \code{colData} and SOM codes in \code{metadata(sce)[[somCodesName]]}.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#'
#' @return A \code{data.frame} of per-cluster statistics.
#'
#' @keywords internal
.buildSOMStats <- function(sce,
    somCodesName = "SOM_codes",
    assay = "exprs") {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]
    colData <- SingleCellExperiment::colData(sce)
    exprs <- SummarizedExperiment::assay(sce, assay)
    stats <- data.frame(row.names = seq_len(nrow(somCodes)))

    for (lv in seq_len(nrow(somCodes))) {
        cells <- exprs[colnames(somCodes), colData$cluster_id == lv,
            drop = FALSE]
        if (ncol(cells) == 0) {
            stats[lv, c("median", "mean", "rdQu", "max", "n", "id")] <-
                c(NA, NA, NA, NA, 0, lv)
            next
        }
        d <- sqrt(colSums((t(cells) - somCodes[lv, , drop = TRUE])^2))
        s1 <- summary(d)
        stats[lv, "median"] <- s1[["Median"]]
        stats[lv, "mean"] <- s1[["Mean"]]
        stats[lv, "rdQu"] <- s1[["3rd Qu."]]
        stats[lv, "max"] <- s1[["Max."]]
        stats[lv, "n"] <- ncol(cells)
        stats[lv, "id"] <- lv
    }
    stats
}


#' Build SOM Cluster Stats
#'
#' Computes per-cluster summary statistics for each marker relative to the SOM
#' code, returning a long-format data frame suitable for downstream plotting.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id} in
#'   \code{colData} and SOM codes in \code{metadata(sce)[[somCodesName]]}.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A long \code{data.frame} with columns \code{cluster}, \code{stat},
#'   \code{variable}, and \code{value}.
#'
#' @keywords internal
.buildSOMClusterStats <- function(sce, somCodesName = "SOM_codes",
    assay = "exprs", verbose = TRUE) {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]
    markers <- colnames(somCodes)
    exprs <- SummarizedExperiment::assay(sce, assay)[markers, , drop = FALSE]
    clusterId <- as.integer(SingleCellExperiment::colData(sce)$cluster_id)
    # Precompute cluster membership as a list of cell indices once.
    idxByCluster <- split(seq_along(clusterId), clusterId)
    result <- lapply(seq_len(nrow(somCodes)), function(lv) {
        if (verbose && lv %% 100 == 0) {
            message("cluster ", lv, " ", format(Sys.time(), "%a %b %d %X %Y"))
        }
        idx <- idxByCluster[[as.character(lv)]]
        cells <- if (length(idx) > 0L) {
            exprs[, idx, drop = FALSE] - somCodes[lv, , drop = TRUE]
        }
        .clusterMarkerSummary(cells, lv, markers)
    })
    do.call(rbind, result)
}

#' Six-Number Summary per Marker of One Cluster
#'
#' @param cells Marker x cell matrix of the cluster's cells minus the node
#'   code, or \code{NULL} for an empty cluster.
#' @param lv Cluster (node) id.
#' @param markers Marker names (rows of \code{cells}).
#'
#' @return A long data frame with \code{cluster}, \code{stat},
#'   \code{variable} and \code{value} (zeros for an empty cluster).
#'
#' @keywords internal
.clusterMarkerSummary <- function(cells, lv, markers) {
    statsNames <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    value <- 0
    if (!is.null(cells)) {
        # matrixStats is much faster than apply() + summary().
        qs <- matrixStats::rowQuantiles(cells,
            probs = c(0, 0.25, 0.5, 0.75, 1), type = 7, drop = FALSE)
        value <- c(qs[, 1], qs[, 2], qs[, 3], matrixStats::rowMeans2(cells),
            qs[, 4], qs[, 5])
        value[is.na(value)] <- 0
    }
    data.frame(
        cluster = rep(lv, length(markers) * length(statsNames)),
        stat = rep(statsNames, each = length(markers)),
        variable = rep(markers, length(statsNames)),
        value = value
    )
}


#' Build SOM Raster Data
#'
#' Builds a data frame representing the SOM grid, including the SOM codes and
#' per-node cell counts.
#'
#' @param sce A \code{SingleCellExperiment} with SOM codes in
#'   \code{metadata(sce)[[somCodesName]]} and optional \code{SOM_stats}.
#'
#' @return A \code{data.frame}.
#'
#' @keywords internal
.buildSOMRasterData <- function(sce, somCodesName = "SOM_codes") {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]
    nDim <- sqrt(nrow(somCodes))
    if (nDim != as.integer(nDim)) {
        stop("SOM grid is not square: nrow(somCodes) is not a perfect square")
    }
    nDim <- as.integer(nDim)
    somRasterData <- data.frame(
        x = rep(seq_len(nDim), nDim),
        y = rep(seq_len(nDim), each = nDim)
    )
    somRasterData <- cbind(somRasterData, somCodes)

    somStats <- S4Vectors::metadata(sce)$SOM_stats
    if (!is.null(somStats) && is.data.frame(somStats) &&
        "n" %in% colnames(somStats)) {
        n_values <- somStats$n
    } else {
        n_values <- rep(NA_integer_, nrow(somCodes))
    }
    somRasterData$N <- n_values
    colnames(somRasterData) <- make.names(colnames(somRasterData))
    somRasterData
}


#' Build Dendrogram Table
#'
#' Creates a table describing the hierarchical clustering dendrogram of the SOM
#' codes, optionally building a collapsible tree network widget.
#'
#' @param sce A \code{SingleCellExperiment} with SOM codes in
#'   \code{metadata(sce)[[somCodesName]]}.
#' @param useColTree Logical indicating whether to build the collapsible tree.
#'
#' @return A list with \code{dend}, \code{dendTable}, and \code{colTree}.
#'
#' @keywords internal
.buildDendTable <- function(sce, somCodesName = "SOM_codes",
    useColTree = FALSE) {
    hc <- stats::hclust(stats::dist(S4Vectors::metadata(sce)[[somCodesName]]))
    dend <- stats::as.dendrogram(hc)
    dendTable <- data.frame(
        parent = NA,
        child = "rt",
        nleaf = .nDendrogramLeaves(dend),
        indexString = ""
    )

    colTree <- NULL
    if (useColTree) {
        # Depth-first walk returning the rows of every sub-tree, in the same
        # order the rows were previously appended.
        recDendTable <- function(dend, parent = "rt", myindexString = "") {
            if (.isDendrogramLeaf(dend)) {
                return(NULL)
            }
            rows <- lapply(seq_along(dend), function(di) {
                indexString <- paste0(myindexString, "[[", di, "]]")
                child <- paste0(parent, ".", di)
                rbind(
                    data.frame(
                        parent = parent,
                        child = child,
                        nleaf = .nDendrogramLeaves(dend[[di]]),
                        indexString = indexString
                    ),
                    recDendTable(dend[[di]], child, indexString)
                )
            })
            do.call(rbind, rows)
        }
        dendTable <- rbind(dendTable, recDendTable(dend))
        colTree <- collapsibleTree::collapsibleTreeNetwork(dendTable,
            attribute = "nleaf", inputId = "node")
    }

    list(dend = dend, dendTable = dendTable, colTree = colTree)
}

#' Test if an Object is a Dendrogram Leaf
#'
#' @param x A dendrogram object.
#'
#' @return Logical.
#'
#' @keywords internal
.isDendrogramLeaf <- function(x) {
    !is.null(attr(x, "leaf")) && isTRUE(attr(x, "leaf"))
}

#' Count Leaves in a Dendrogram
#'
#' @param x A dendrogram object.
#'
#' @return Integer.
#'
#' @keywords internal
.nDendrogramLeaves <- function(x) {
    if (.isDendrogramLeaf(x)) {
        return(1L)
    }
    sum(vapply(x, .nDendrogramLeaves, integer(1)))
}


#' Prepare Cluster Selector Outputs
#'
#' Wrapper that takes a post-SOM \code{SingleCellExperiment}, builds all
#' supporting data structures required by \code{\link{clusterSelector}}, and
#' saves a \code{sce.shiny.RData} bundle compatible with the original
#' \code{prepCytoApp} output.
#'
#' @param currentRoot Character string specifying the project root directory.
#'   Output files are written to \code{currentRoot/_data/}.
#' @param sce A \code{SingleCellExperiment} with SOM clustering. If
#'   \code{NULL}, the function reads \code{currentRoot/_data/sce.RDS}.
#' @param useColTree Logical indicating whether to create a collapsible tree.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param clusterColName Column name containing cluster IDs.
#' @param total_cells_to_sample Number of cells to subsample for the Shiny app.
#' @param dList Optional list of marker pairs for 2D plots.
#' @param assay Name of the assay to use.
#' @param verbose Logical indicating whether to print progress messages.
#' @param force Logical indicating whether to force recomputation.
#' @param outputList Optional output list to include in the saved bundle.
#' @param seed Random seed for subsampling.
#'
#' @return The path \code{currentRoot}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsample the SCE using \code{prepClusterSelectorData()}, or reuse a
#'     previously subsampled SCE file if \code{subsampledSCEFile} is provided.
#'   \item Update \code{cluster_id} to a factor with all SOM clusters as levels.
#'   \item Build the cluster-by-patient table.
#'   \item Recompute missing SOM code means.
#'   \item Build SOM summary stats and cluster-level stats.
#'   \item Build SOM raster data.
#'   \item Build hierarchical clustering dendrogram and table.
#'   \item Save \code{sce.shiny.RData}.
#' }
#'
#' @param subsampledSCEFile Optional path to a previously subsampled SCE object
#'   (for example \code{sce.subsampled.RDS}). If provided and the file exists,
#'   it is reused instead of subsampling the full SCE. This avoids the memory
#'   cost of subsampling large datasets.
#'
#' @examples
#' root <- file.path(tempdir(), "cysa-example")
#' dir.create(file.path(root, "_data"), recursive = TRUE, showWarnings = FALSE)
#' # A square SOM grid (here 4 x 4 nodes) is required for the raster data.
#' sce <- CySA_example_sce(n_cells = 200, n_nodes = 16)
#' saveRDS(sce, file.path(root, "_data", "sce.RDS"))
#' prepClusterSelector(root, sce = sce, total_cells_to_sample = 100,
#'     verbose = FALSE)
#' list.files(file.path(root, "_data"))
#'
#' @export
prepClusterSelector <- function(currentRoot, sce = NULL, useColTree = FALSE,
    somCodesName = "SOM_codes", clusterColName = "cluster_id",
    total_cells_to_sample = 100000, dList = NULL, assay = "exprs",
    verbose = TRUE, force = FALSE, outputList = NULL, seed = 123,
    subsampledSCEFile = NULL) {
    say <- function(...) if (verbose) message(...)
    somFile <- file.path(currentRoot, "_data", "sce.RDS")
    outFile <- file.path(currentRoot, "_data", "sce.shiny.RData")
    if (!file.exists(somFile)) {
        stop("SCE file not found: ", somFile)
    }
    if (file.exists(outFile) && !force) {
        say("Output exists and force = FALSE; returning currentRoot")
        return(currentRoot)
    }
    if (is.null(sce)) {
        say("Reading SCE")
        sce <- readRDS(somFile)
    }
    sub <- .prepSubsample(sce, subsampledSCEFile, somFile, dList,
        total_cells_to_sample, somCodesName, assay, seed, say)
    sce_subsampled <- sub$sce_subsampled
    dList <- sub$dList
    say("Building cluster patient table")
    sce$cluster_id <- factor(
        as.character(SingleCellExperiment::colData(sce)[[clusterColName]]),
        levels = seq_len(nrow(S4Vectors::metadata(sce)[[somCodesName]]))
    )
    clusterPatientTable <- table(as.data.frame(
        SingleCellExperiment::colData(sce)[, c("sample_id", clusterColName)]
    ))
    sce <- .prepSomSummaries(sce, somCodesName, assay, verbose, say)
    say("Building SOM raster data")
    somRasterData <- .buildSOMRasterData(sce, somCodesName = somCodesName)
    somRasterObj <- NULL
    say("Building dendrogram table")
    dendRes <- .buildDendTable(sce, somCodesName = somCodesName,
        useColTree = useColTree)
    # save() stores objects under their variable names: the bundle holds
    # dend_obj and dend_table, which is what the loaders expect.
    dend_obj <- dendRes$dend
    dend_table <- dendRes$dendTable
    say("Saving outputs")
    save(sce, sce_subsampled, outputList, dList, dend_obj, dend_table,
        clusterPatientTable, somRasterData, somRasterObj,
        file = outFile, compress = FALSE)
    say("Saved: ", outFile)
    currentRoot
}

#' Subsampled SCE for prepClusterSelector()
#'
#' Reuses \code{subsampledSCEFile} when it exists (restoring its assay if it
#' was saved without one), otherwise subsamples with
#' \code{prepClusterSelectorData()}.
#'
#' @inheritParams prepClusterSelector
#' @param somFile Path of the saved full SCE (cache key of the subsample).
#' @param say Progress message function.
#'
#' @return A list with \code{sce_subsampled} and \code{dList}.
#'
#' @keywords internal
.prepSubsample <- function(sce, subsampledSCEFile, somFile, dList,
    total_cells_to_sample, somCodesName, assay, seed, say) {
    if (!is.null(subsampledSCEFile) && file.exists(subsampledSCEFile)) {
        say("Using provided subsampled SCE file: ", subsampledSCEFile)
        if (is.null(dList)) {
            dList <- .defaultDList(sce)
        }
        sub <- readRDS(subsampledSCEFile)
        if (!is.null(sub) &&
            length(SummarizedExperiment::assays(sub)) == 0) {
            say("Restoring assay data into subsampled SCE")
        }
        return(list(sce_subsampled = .restoreSubsampledAssay(sub, sce, assay),
            dList = dList))
    }
    say("Subsampling SCE")
    prep <- prepClusterSelectorData(
        sce = sce, somFile = somFile, dList = dList,
        total_cells_to_sample = total_cells_to_sample,
        somCodesName = somCodesName, assay = assay, seed = seed
    )
    list(sce_subsampled = prep$sce_subsampled, dList = prep$dList)
}

#' SOM Codes, Node Statistics and Cluster Statistics
#'
#' @inheritParams prepClusterSelector
#' @param say Progress message function.
#'
#' @return \code{sce} with refreshed SOM codes, \code{SOM_stats} (built only
#'   when missing) and \code{som_clster_stats} in its metadata.
#'
#' @keywords internal
.prepSomSummaries <- function(sce, somCodesName, assay, verbose, say) {
    say("Updating SOM codes")
    sce <- .buildSOMCodes(sce, somCodesName = somCodesName, assay = assay)
    if (is.null(S4Vectors::metadata(sce)$SOM_stats)) {
        say("Building SOM stats")
        S4Vectors::metadata(sce)$SOM_stats <- .buildSOMStats(sce,
            somCodesName = somCodesName, assay = assay)
    } else {
        say("Reusing existing SOM_stats")
    }
    say("Building SOM cluster stats")
    S4Vectors::metadata(sce)[["som_clster_stats"]] <- .buildSOMClusterStats(
        sce, somCodesName = somCodesName, assay = assay, verbose = verbose
    )
    sce
}
