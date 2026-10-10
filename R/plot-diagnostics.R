# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# plot-diagnostics.R ----
# Pure plotting/data functions for the sample-level QC diagnostics (sample
# MDS, per-sample marker heatmap, marker non-redundancy scores) and the
# dendrogram-ordered cluster heatmap added to CySA on top of
# vignettes/catalyst_workflow.Rmd, which reproduces the equivalent
# CATALYST::pbMDS()/plotExprHeatmap()/plotNRS()/plotMultiHeatmap() figures
# with the reference implementation. No Shiny reactivity lives here so these
# can be tested and reused independently of the app.

#' Subset a SingleCellExperiment to Cells in the Given Saved Groups
#'
#' Lets the diagnostics/ordered-heatmap panels restrict their computation to
#' a user-chosen selection of saved clusters/metaclusters (CySA's
#' \code{rv$outputList}) instead of always using every cell.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{cluster_id} colData
#'   column (SOM node ids).
#' @param outputList Named list of SOM node ids per saved group.
#' @param groups Names of \code{outputList} to keep; \code{NULL}, empty, or
#'   every group selected keeps the whole dataset unfiltered.
#'
#' @return \code{sce}, subset to cells whose node belongs to one of
#'   \code{groups}.
#'
#' @keywords internal
.subsetSceByGroups <- function(sce, outputList, groups = NULL) {
    if (is.null(groups) || length(groups) == 0L ||
        setequal(groups, names(outputList))) {
        return(sce)
    }
    nodeIds <- unlist(lapply(outputList[groups], as.character),
        use.names = FALSE)
    keep <- as.character(SingleCellExperiment::colData(sce)$cluster_id) %in%
        nodeIds
    sce[, keep]
}

#' Per-Sample Pseudobulk Marker Medians
#'
#' @param sce A \code{SingleCellExperiment}; the first assay is used.
#' @param markers Marker (row) names to include.
#'
#' @return A samples x markers matrix of per-sample median marker expression.
#'
#' @keywords internal
.samplePseudobulk <- function(sce, markers) {
    exprMat <- SummarizedExperiment::assay(sce, 1)[markers, , drop = FALSE]
    sampleId <- SingleCellExperiment::colData(sce)$sample_id
    samples <- unique(sampleId)
    pb <- vapply(samples, function(s) {
        matrixStats::rowMedians(exprMat[, sampleId == s, drop = FALSE])
    }, numeric(length(markers)))
    pb <- t(pb)
    rownames(pb) <- as.character(samples)
    colnames(pb) <- markers
    pb
}

#' Attach an experiment_info Column for Colouring
#'
#' @param df Data frame with a \code{sample_id} column.
#' @param sce A \code{SingleCellExperiment} whose
#'   \code{metadata()$experiment_info} supplies the column.
#' @param color_by Column name in \code{experiment_info}, or \code{NULL}/
#'   \code{"none"}.
#'
#' @return \code{df}, with a \code{color_by} column added when possible.
#'
#' @keywords internal
.attachColorBy <- function(df, sce, color_by) {
    if (is.null(color_by) || identical(color_by, "none")) {
        return(df)
    }
    info <- S4Vectors::metadata(sce)$experiment_info
    if (is.null(info) || !color_by %in% colnames(info)) {
        return(df)
    }
    df$color_by <- info[match(df$sample_id, info$sample_id), color_by]
    df
}

#' Pseudobulk MDS Plot of the Samples
#'
#' Equivalent of \code{CATALYST::pbMDS()}: multidimensional scaling of the
#' per-sample pseudobulk marker medians.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param markers Marker names used to compute the pseudobulk distances.
#' @param color_by Optional \code{experiment_info} column used to colour
#'   samples.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plotSampleMDS <- function(sce, markers, color_by = NULL) {
    pb <- .samplePseudobulk(sce, markers)
    if (nrow(pb) < 3L) {
        stop("plotSampleMDS() needs at least 3 samples.")
    }
    xy <- stats::cmdscale(stats::dist(pb), k = 2)
    df <- data.frame(
        sample_id = rownames(pb), mds1 = xy[, 1], mds2 = xy[, 2]
    )
    df <- .attachColorBy(df, sce, color_by)
    p <- ggplot(df, aes(x = mds1, y = mds2, label = sample_id))
    p <- if ("color_by" %in% colnames(df)) {
        p + geom_point(aes(color = color_by), size = 3)
    } else {
        p + geom_point(size = 3)
    }
    p + geom_text(vjust = -1.2, size = 3) +
        labs(x = "MDS 1", y = "MDS 2", color = color_by %||% NULL) +
        theme_minimal()
}

#' Per-Sample Marker Expression Heatmap
#'
#' Equivalent of \code{CATALYST::plotExprHeatmap()}: median marker expression
#' per sample, with samples and markers each ordered by their own
#' hierarchical clustering.
#'
#' @inheritParams plotSampleMDS
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plotSampleExprHeatmap <- function(sce, markers) {
    pb <- .samplePseudobulk(sce, markers)
    sampleOrd <- if (nrow(pb) >= 2L) {
        rownames(pb)[stats::hclust(stats::dist(pb))$order]
    } else {
        rownames(pb)
    }
    markerOrd <- if (ncol(pb) >= 2L) {
        colnames(pb)[stats::hclust(stats::dist(t(pb)))$order]
    } else {
        colnames(pb)
    }
    df <- as.data.frame(pb)
    df$sample_id <- rownames(pb)
    long <- tidyr::pivot_longer(df, cols = tidyselect::all_of(markers),
        names_to = "marker", values_to = "value")
    long$sample_id <- factor(long$sample_id, levels = sampleOrd)
    long$marker <- factor(long$marker, levels = markerOrd)
    ggplot(long, aes(x = marker, y = sample_id, fill = value)) +
        geom_tile() +
        scale_fill_viridis() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Marker Non-Redundancy Scores
#'
#' Per-sample marker non-redundancy score (NRS), the published Nowicka
#' \emph{et al.} (2017) definition: for each sample, the squared loading of
#' each marker on the first three principal components, weighted by the
#' percentage of variance each component explains.
#'
#' @inheritParams plotSampleMDS
#'
#' @return A data frame with \code{sample_id}, \code{marker} and \code{nrs}.
#'
#' @keywords internal
computeNRS <- function(sce, markers) {
    exprMat <- t(SummarizedExperiment::assay(sce, 1)[markers, , drop = FALSE])
    sampleId <- SingleCellExperiment::colData(sce)$sample_id
    samples <- unique(sampleId)
    rows <- lapply(samples, function(s) {
        m <- exprMat[sampleId == s, , drop = FALSE]
        if (nrow(m) < 3L) {
            return(NULL)
        }
        pr <- stats::prcomp(m, center = TRUE, scale. = FALSE)
        nPc <- min(3L, ncol(pr$rotation))
        pctVar <- (pr$sdev[seq_len(nPc)]^2) / sum(pr$sdev^2)
        nrs <- as.numeric(pr$rotation[, seq_len(nPc), drop = FALSE]^2 %*%
            pctVar)
        data.frame(sample_id = s, marker = markers, nrs = nrs)
    })
    do.call(rbind, rows)
}

#' Marker Non-Redundancy Score Plot
#'
#' Equivalent of \code{CATALYST::plotNRS()}.
#'
#' @inheritParams plotSampleMDS
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plotNRS <- function(sce, markers, color_by = NULL) {
    nrsDf <- computeNRS(sce, markers)
    ordStat <- stats::aggregate(nrs ~ marker, data = nrsDf,
        FUN = stats::median)
    ord <- ordStat$marker[order(-ordStat$nrs)]
    nrsDf$marker <- factor(nrsDf$marker, levels = ord)
    nrsDf <- .attachColorBy(nrsDf, sce, color_by)
    p <- ggplot(nrsDf, aes(x = marker, y = nrs))
    p <- if ("color_by" %in% colnames(nrsDf)) {
        p + geom_jitter(aes(color = color_by), width = 0.2, size = 2)
    } else {
        p + geom_jitter(width = 0.2, size = 2)
    }
    p + geom_boxplot(alpha = 0.3, outlier.shape = NA) +
        labs(x = NULL, y = "Non-redundancy score",
            color = color_by %||% NULL) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Dendrogram-Ordered Cluster Heatmap
#'
#' Equivalent of \code{CATALYST::plotMultiHeatmap(bars = TRUE)}: marker
#' medians per cluster, ordered by the existing SOM dendrogram's leaf order,
#' with a cell-count bar alongside.
#'
#' @param somCodes Nodes x markers matrix (e.g. \code{metadata(sce)$SOM_codes}).
#' @param somStats Data frame with node \code{id} and cell count \code{n}
#'   (e.g. \code{metadata(sce)$SOM_stats}), or \code{NULL} to omit the bar.
#' @param dend Dendrogram object (leaf order used to reorder clusters).
#' @param dendTable Optional data frame with \code{id}/\code{label}, mapping
#'   dendrogram leaf labels back to cluster ids.
#' @param nodeIds Cluster/SOM-node ids to display; \code{NULL} shows every
#'   cluster in \code{somCodes}.
#' @param markers Marker (column) names to display; \code{NULL} shows every
#'   marker in \code{somCodes}.
#'
#' @return A \code{ggplot} object, or a \code{cowplot} grid when
#'   \code{somStats} is supplied.
#'
#' @keywords internal
plotOrderedClusterHeatmap <- function(somCodes, somStats, dend,
    dendTable = NULL, nodeIds = NULL, markers = NULL) {
    g <- dendextend::as.ggdend(dend)
    ord <- as.character(g$labels$label)
    if (!is.null(dendTable) &&
        all(c("id", "label") %in% colnames(dendTable))) {
        ord <- as.character(dendTable$id[match(ord, dendTable$label)])
    }
    ord <- intersect(ord, rownames(somCodes))
    if (!is.null(nodeIds)) {
        ord <- intersect(ord, as.character(nodeIds))
    }
    if (length(ord) == 0L) {
        stop("No clusters left to show; select at least one cluster/group.")
    }
    cols <- if (is.null(markers)) {
        colnames(somCodes)
    } else {
        intersect(markers, colnames(somCodes))
    }
    if (length(cols) == 0L) {
        stop("No markers left to show; select at least one marker.")
    }
    somCodes <- somCodes[ord, cols, drop = FALSE]
    df <- as.data.frame(somCodes)
    df$cluster_id <- factor(ord, levels = ord)
    long <- tidyr::pivot_longer(df, cols = -"cluster_id",
        names_to = "marker", values_to = "value")
    heat <- ggplot(long, aes(x = marker, y = cluster_id, fill = value)) +
        geom_tile() +
        scale_fill_viridis() +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6)
        )
    if (is.null(somStats) || !all(c("id", "n") %in% colnames(somStats))) {
        return(heat)
    }
    nDf <- data.frame(
        cluster_id = factor(ord, levels = ord),
        n = somStats$n[match(as.integer(ord), somStats$id)]
    )
    bar <- ggplot(nDf, aes(x = n, y = cluster_id)) +
        geom_col(fill = "grey40") +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(), axis.title.y = element_blank()
        ) +
        labs(x = "n cells")
    cowplot::plot_grid(heat, bar, ncol = 2, rel_widths = c(3, 1), align = "h")
}
