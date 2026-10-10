# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# functions used in shiny app
#' @importFrom tidyselect all_of
# safe_event_data ----
# Wrapper around plotly::event_data that suppresses the "not registered"
# warnings unless verbose mode is enabled.
safe_event_data <- function(event, source = "all", verbose = FALSE) {
    if (verbose) {
        return(plotly::event_data(event, source = source))
    }
    withCallingHandlers(
        plotly::event_data(event, source = source),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("is not registered", msg) ||
                grepl("event_register", msg)) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

# Alias with dot prefix used by extracted server helpers.
.safeEventData <- safe_event_data

# quiet_ggplotly ----
# Wrapper around ggplotly that suppresses the harmless ggplot2 warning about
# the customdata aesthetic, which is consumed by plotly after conversion.
quiet_ggplotly <- function(p, ...) {
    result <- withCallingHandlers(
        ggplotly(p, ...),
        warning = function(w) {
            if (.isHarmlessPlotlyWarning(conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        },
        message = function(m) invokeRestart("muffleMessage")
    )
    # Remove ggplot2-internal 'colour' attribute from every trace.
    result$x$data <- lapply(result$x$data, function(tr) {
        tr[["colour"]] <- NULL
        tr
    })
    result <- result %>%
        plotly::layout(autosize = TRUE) %>%
        plotly::config(responsive = TRUE)
    # Force early evaluation of the lazy plotly object. plotly_build() is what
    # actually resolves trace type/mode and emits the "No trace type
    # specified"/"No scatter mode specified" messages -- muffling them here
    # means Shiny's later renderPlotly()/toJSON step (or a test accessing
    # output$xxx) just serializes an already-resolved object and has nothing
    # left to emit.
    withCallingHandlers(
        plotly::plotly_build(result),
        message = function(m) invokeRestart("muffleMessage")
    )
}

#' Known Harmless Warnings of the ggplot to plotly Conversion
#'
#' @param msg Warning message.
#'
#' @return \code{TRUE} for the customdata/key aesthetic, aspect-ratio,
#'   not-yet-implemented and unknown-attribute warnings.
#'
#' @keywords internal
.isHarmlessPlotlyWarning <- function(msg) {
    patterns <- c(
        "unknown aesthetics.*(customdata|key)",
        "Aspect ratios aren't yet implemented",
        "has yet to be implemented in plotly",
        "don't have these attributes"
    )
    any(vapply(patterns, grepl, logical(1), x = msg, ignore.case = TRUE))
}

# highlight_df = function ----
highlight_df <- function(x, y, rs, somCodesName = "SOM_codes", metaD = NULL) {
    if (is.null(metaD)) stop("'metaD' must be provided")
    if (is.null(rs) || length(rs) == 0L) {
        return(data.frame(x = numeric(0), y = numeric(0), id = integer(0)))
    }

    # Filter out NA and out-of-bounds indices
    somCodes <- metaD[[somCodesName]]
    rs_valid <- rs[!is.na(rs) & rs > 0L & rs <= nrow(somCodes)]

    if (length(rs_valid) == 0L) {
        return(data.frame(x = numeric(0), y = numeric(0), id = integer(0)))
    }

    data.frame(
        x = somCodes[rs_valid, x],
        y = somCodes[rs_valid, y],
        id = rs_valid
    )
}


countBarPlotFunc <- function(rs, clusterPatientTable, cst, sce, outputList,
    groupsInput) {
    rs <- as.integer(intersect(rs, colnames(clusterPatientTable)))
    rSums <- data.frame(
        counts = rowSums(clusterPatientTable[, rs, drop = FALSE])
    )
    rSums$id <- rownames(clusterPatientTable)
    info <- S4Vectors::metadata(sce)$experiment_info
    dfm <- .countComparison(rSums, cst, info, clusterPatientTable, outputList)
    if ("group" %in% colnames(info)) {
        dfm <- merge(dfm, info, by.x = "id", by.y = "sample_id")
        dfm$id <- factor(dfm$id,
            levels = info[order(info$group), "sample_id"])
    }
    dfm <- as_tibble(.tagSampleGroups(dfm, groupsInput))
    ggplot(dfm, aes(x = reorder(id, counts), y = counts)) +
        geom_bar(aes(fill = interaction(variable, groups)),
            stat = "identity", position = "dodge") +
        coord_flip() +
        guides(fill = guide_legend(title = "groups"))
}

#' Selected Counts With an Optional Comparison Series
#'
#' @param rSums Data frame with per-sample \code{counts} and \code{id}.
#' @param cst Comparison: "none", a numeric \code{experiment_info} column or
#'   a saved group name.
#' @param info The \code{experiment_info} data frame.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param outputList Named list of saved groups.
#'
#' @return A long data frame with \code{id}, \code{variable} and
#'   \code{counts}.
#'
#' @keywords internal
.countComparison <- function(rSums, cst, info, clusterPatientTable,
    outputList) {
    if (cst == "none") {
        rSums$variable <- "counts"
        return(rSums)
    }
    if (cst %in% colnames(info)) {
        numCols <- unlist(lapply(info, is.numeric), use.names = FALSE)
        eI <- apply(info[, numCols, drop = FALSE], 2, as.numeric) %>%
            as.data.frame()
        rownames(eI) <- info[, "sample_id"]
        if (any(is.na(eI))) {
            stop("NAs produced when converting experiment_info to numeric")
        }
        comp <- eI[rSums$id, cst]
    } else if (cst %in% names(outputList)) {
        nodes <- .validNodeCols(outputList[[cst]], clusterPatientTable)
        comp <- rowSums(clusterPatientTable[, nodes, drop = FALSE])
    } else {
        stop("Comparison column '", cst,
            "' not found in experiment_info or outputList")
    }
    rSums <- cbind(rSums, comp)
    colnames(rSums) <- c("counts", "id", cst)
    dfm <- reshape2::melt(rSums[, c("id", cst, "counts")], id.vars = 1)
    colnames(dfm) <- c("id", "variable", "counts")
    dfm
}

#' Tag Samples With Their t-Test Group
#'
#' @param dfm Data frame with a sample \code{id} column.
#' @param groupsInput List with \code{group1} and \code{group2} sample ids,
#'   or \code{NULL}.
#'
#' @return \code{dfm} with a factor column \code{groups} ("group1",
#'   "group2" or "other").
#'
#' @keywords internal
.tagSampleGroups <- function(dfm, groupsInput) {
    dfm$groups <- "other"
    if (!is_empty(groupsInput)) {
        if (!is_empty(groupsInput$group1)) {
            dfm$groups[dfm$id %in% groupsInput$group1] <- "group1"
        }
        if (!is_empty(groupsInput$group2)) {
            dfm$groups[dfm$id %in% groupsInput$group2] <- "group2"
        }
    }
    dfm$groups <- as.factor(dfm$groups)
    dfm
}


# Deprecated alias kept for backward compatibility.
# Prefer .computeRelativeCounts() in new code.
compute_relative_counts <- function(clusterPatientTable, rs, relativeToCol,
    expInfo, outputList) {
    .computeRelativeCounts(clusterPatientTable, rs, relativeToCol, expInfo,
        outputList)
}

PercentBarPlotFunc <- function(sce, relativeToCol, clusterPatientTable, rs,
    outputList, group, groupsInput) {
    info <- S4Vectors::metadata(sce)$experiment_info
    expInfo <- info
    rownames(expInfo) <- expInfo$sample_id
    expInfo <- expInfo[, !colnames(expInfo) %in%
        c("sample_nr", "sample_id", "sample"), drop = FALSE]
    rSums <- data.frame(Percent = .computeRelativeCounts(
        clusterPatientTable, rs, relativeToCol, expInfo, outputList
    ))
    rSums$id <- rownames(clusterPatientTable)
    if ("group" %in% colnames(info)) {
        rSums <- merge(rSums, info, by.x = "id", by.y = "sample_id")
        rSums <- rSums %>% group_by(group)
        rSums$id <- factor(rSums$id,
            levels = rSums[order(rSums$group), "id"][[1]])
    }
    rSums <- as_tibble(.tagSampleGroups(rSums, groupsInput))
    ggplot(rSums, aes(x = reorder(id, Percent), y = Percent, fill = groups)) +
        geom_col() +
        coord_flip()
}


ggsomPlot <- function(pp1, plotIdx, rs, dimSelection,
    somCodesName = "SOM_codes", sce, metaD = S4Vectors::metadata(sce),
    xlim = NULL, ylim = NULL) {
    ds <- dimSelection[[plotIdx]]
    newData <- highlight_df(ds$dims[1], ds$dims[2], rs, somCodesName,
        metaD = metaD)
    p3 <- pp1 + geom_point(data = newData, aes(x = `x`, y = `y`),
        color = "red", size = 0.3)
    # coord_cartesian() clips instead of dropping points. Explicit limits
    # win; otherwise the stored channel limits, or the user's zoom.
    num2 <- function(v) as.numeric(v[c(1, 2)])
    p3 + if (!is.null(xlim) && !is.null(ylim)) {
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
    } else if (is.null(ds$xzoom[1])) {
        ggplot2::coord_cartesian(xlim = num2(ds$xlim), ylim = num2(ds$ylim))
    } else {
        ggplot2::coord_cartesian(xlim = num2(ds$xzoom), ylim = num2(ds$yzoom))
    }
}

somPlot <- function(pp1, plotIdx, rs, colorbyGroups, showGroups,
    dimSelection = NULL, somCodesName = "SOM_codes", sce,
    metaD = S4Vectors::metadata(sce), outputList = list(),
    projectionDf = NULL, xlim = NULL, ylim = NULL, source = NULL,
    agreement = NULL) {
    if (is.null(pp1)) {
        return(NULL)
    }
    # The caller (somDataMain) passes a single-element dimSelection.
    ds <- dimSelection[[1L]]
    # The agreement view needs the per-node projection data, so it routes
    # through drawProjection() even when group colouring is off.
    p3 <- if (showGroups || !is.null(agreement)) {
        req(projectionDf)
        drawProjection(projectionDf, rs, colorbyGroups = colorbyGroups,
            sce = sce, outputList = outputList, agreement = agreement)
    } else {
        pp1 + geom_point(
            data = highlight_df(ds$dims[1], ds$dims[2], rs, somCodesName,
                metaD = metaD),
            aes(x = x, y = y), color = "red", size = 0.3
        )
    }
    if (is.null(p3)) {
        return(NULL)
    }
    # Limits are applied to the ggplot before the plotly conversion.
    p3 <- p3 + .somCoord(ds, xlim, ylim)
    # Use an explicit source when supplied; otherwise derive it from the
    # plot index.
    if (plotIdx > 5L) plotIdx <- 6L
    plotSource <- if (!is.null(source)) source else paste0("somData", plotIdx)
    quiet_ggplotly(p3, source = plotSource, tooltip = "") %>%
        layout(showlegend = FALSE, dragmode = "select") %>%
        event_register("plotly_selected") %>%
        event_register("plotly_relayout")
}

#' Axis Limits of a SOM Plot
#'
#' The user's zoom wins; otherwise the supplied \code{xlim}/\code{ylim};
#' otherwise the channel limits stored in the dimension selection.
#'
#' @param ds One element of the dimension selection (\code{dims},
#'   \code{xlim}, \code{ylim}, \code{xzoom}, \code{yzoom}).
#' @param xlim,ylim Optional explicit limits.
#'
#' @return A \code{coord_cartesian}, or \code{NULL} when no limits apply.
#'
#' @keywords internal
.somCoord <- function(ds, xlim = NULL, ylim = NULL) {
    num2 <- function(v) as.numeric(v[c(1, 2)])
    if (!is.null(ds$xzoom[1])) {
        return(ggplot2::coord_cartesian(xlim = num2(ds$xzoom),
            ylim = num2(ds$yzoom)))
    }
    if (!is.null(xlim) && !is.null(ylim)) {
        return(ggplot2::coord_cartesian(xlim = xlim, ylim = ylim))
    }
    if (!is.null(ds$xlim) && !is.null(ds$ylim)) {
        return(ggplot2::coord_cartesian(xlim = num2(ds$xlim),
            ylim = num2(ds$ylim)))
    }
    NULL
}


tsneFunc <- function(dimRedSelection, perplexity, sce,
    somCodesName = "SOM_codes") {
    dimRedCols <- dimRedSelection
    mat <- S4Vectors::metadata(sce)[[somCodesName]][, dimRedCols, drop = FALSE]
    max_perplexity <- max(1, floor((nrow(mat) - 1) / 3))
    perplexity <- min(as.integer(perplexity), max_perplexity)
    Rtsne::Rtsne(mat, perplexity = perplexity)
}


plotViolinFunc <- function(sce, somCodesName = "SOM_codes", upsetSelection,
    outputList, violinSelection) {
    if (length(upsetSelection) < 3) {
        upsetSelection <- names(outputList)
    }
    data <- .violinData(sce@metadata[[somCodesName]], upsetSelection,
        outputList, violinSelection)
    if (nrow(data) > 0L && nrow(data) <= 10L) {
        message("please check that all upsetSelection are in outputList: ",
            "no violin data to plot")
        return(NULL)
    }
    .violinGgplot(data, "No groups with \u2265 2 SOM nodes selected")
}

#' Long-Format SOM Code Data of Saved Groups for Violin Plots
#'
#' @param somCodes SOM codes matrix (nodes x markers).
#' @param groups Names of the groups in \code{outputList} to include.
#' @param outputList Named list of SOM node ids per group.
#' @param markers Markers to include (intersected with the SOM columns).
#'
#' @return A data frame with \code{somNode}, \code{marker}, \code{expr} and
#'   \code{grpName} (factors in marker / group order); groups with fewer
#'   than two nodes are skipped.
#'
#' @keywords internal
.violinData <- function(somCodes, groups, outputList, markers) {
    markers <- intersect(colnames(somCodes), markers)
    nNodes <- nrow(somCodes)
    parts <- lapply(groups, function(na) {
        rows <- outputList[[na]]
        # Drop the placeholder 0 and any id a stale saved group carries
        # beyond the current SOM's node range.
        rows <- rows[!is.na(rows) & rows >= 1L & rows <= nNodes]
        if (length(rows) < 2) {
            return(NULL)
        }
        wide <- as.data.frame(somCodes[rows, markers, drop = FALSE])
        wide$somNode <- factor(rows, levels = seq_len(nNodes))
        long <- tidyr::pivot_longer(wide,
            cols = tidyselect::all_of(markers),
            names_to = "marker", values_to = "expr"
        )
        long$grpName <- na
        long
    })
    data <- as.data.frame(data.table::rbindlist(parts))
    if (nrow(data) > 0L) {
        data$marker <- factor(data$marker, levels = markers)
        data$grpName <- factor(data$grpName, levels = groups)
    }
    data
}

#' Violin Plot of SOM Codes per Marker, Faceted by Group
#'
#' @param data Output of \code{.violinData()}.
#' @param emptyTitle Title of the placeholder plot when \code{data} has no
#'   rows.
#' @param xVar \code{"marker"}: markers on the x-axis, one facet row per
#'   group; \code{"grpName"}: groups on the x-axis, one row per marker.
#'
#' @return A ggplot.
#'
#' @keywords internal
.violinGgplot <- function(data, emptyTitle, xVar = c("marker", "grpName")) {
    xVar <- match.arg(xVar)
    if (nrow(data) == 0L) {
        return(ggplot() + theme_void() + labs(title = emptyTitle))
    }
    byMarker <- xVar == "marker"
    facetVar <- if (byMarker) "grpName" else "marker"
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(
        length(unique(data[[xVar]]))
    )
    ggplot(data, aes(factor(.data[[xVar]]), .data$expr,
        fill = .data[[xVar]])) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_y_continuous(expand = c(0, 0), position = "right",
            labels = function(x) {
                c(rep(x = "", times = length(x) - 2), x[length(x) - 1], "")
            }) +
        facet_grid(rows = vars(.data[[facetVar]]), scales = "free",
            switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(
            legend.position = "none", panel.spacing = unit(0, "lines"),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = NA, color = "black"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.y.left = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) +
        scale_fill_manual(values = mycolors) +
        ggtitle(if (byMarker) "Marker on x-axis" else "grp name on x-axis") +
        xlab(if (byMarker) "Marker" else "Groups") +
        ylab("Expression Level")
}



#' Build a Projection Data Frame for drawProjection
#'
#' Constructs the per-SOM-node data frame that \code{drawProjection} expects:
#' the two requested channel columns come first, followed by all SOM_stats
#' columns.  Reads coordinates directly from
#' \code{metadata(sce)[[somCodesName]]}
#' instead of parsing \code{ggplot_build()} output, which is both faster and
#' independent of which aesthetics the base plot happens to expose.
#'
#' @param pp1 Base ggplot object (kept in signature for call-site compatibility;
#'   no longer used internally).
#' @param plotIdx Index into \code{dimSelection}.
#' @param dimSelection List of per-plot dimension specs.
#' @param sce Full \code{SingleCellExperiment}.
#' @param somCodesName Name of the SOM codes metadata slot
#'   (default \code{"SOM_codes"}).
#'
#' @return A data frame or \code{NULL} on failure.
#'
#' @keywords internal
buildProjectionDf <- function(pp1, plotIdx, dimSelection, sce,
    somCodesName = "SOM_codes") {
    shiny::req(pp1, dimSelection)
    md <- S4Vectors::metadata(sce)
    somCodes <- md[[somCodesName]]
    somStats <- md[["SOM_stats"]]
    if (is.null(somCodes)) {
        warning("buildProjectionDf: '", somCodesName,
            "' not found in metadata(sce).")
        return(NULL)
    }
    if (is.null(somStats)) {
        warning("buildProjectionDf: 'SOM_stats' not found in metadata(sce).")
        return(NULL)
    }
    if (!"id" %in% colnames(somStats)) {
        somStats$id <- seq_len(nrow(somStats))
    }
    chNames <- .resolveProjectionChannels(dimSelection[[plotIdx]]$dims,
        colnames(somCodes))
    # Row i of the SOM codes is SOM node i.
    coordDf <- data.frame(id = seq_len(nrow(somCodes)),
        ch1 = somCodes[, chNames[1L]], ch2 = somCodes[, chNames[2L]],
        stringsAsFactors = FALSE, check.names = FALSE)
    colnames(coordDf)[2:3] <- chNames
    df2 <- dplyr::left_join(somStats, coordDf, by = "id")
    df2 <- df2[order(df2$id), ]
    df2[is.na(df2)] <- 0
    # Channel columns first: drawProjection() reads names(df)[1:2].
    coordIdx <- match(chNames, names(df2))
    df2[, c(coordIdx, setdiff(seq_along(df2), coordIdx)), drop = FALSE]
}

#' Replace Unknown Projection Channels
#'
#' @param chNames Requested channel pair.
#' @param available Columns of the SOM codes.
#'
#' @return \code{chNames}, unknown channels replaced (with a warning) by the
#'   first available ones.
#'
#' @keywords internal
.resolveProjectionChannels <- function(chNames, available) {
    missingCh <- setdiff(chNames, available)
    if (length(missingCh) > 0L) {
        warning(
            "buildProjectionDf: channel(s) not in SOM codes - falling back: ",
            paste(missingCh, collapse = ", ")
        )
        chNames[!chNames %in% available] <-
            head(available, 2L)[seq_along(missingCh)]
    }
    chNames
}


#' kNN Agreement Between a 2D View and the Full SOM Codes
#'
#' For every SOM node, compares its \code{k} nearest neighbours in the
#' displayed 2D marker pair with its \code{k} nearest neighbours in the
#' full SOM-code space. The score is the Jaccard index of the two
#' neighbour sets: 1 means the 2D view preserves the node's local
#' neighbourhood, values near 0 mean the pair distorts it.
#'
#' @param codes Numeric matrix of SOM codes (nodes x markers).
#' @param xy Numeric matrix of 2D coordinates, one row per SOM node, in
#'   the same node order as \code{codes}.
#' @param k Number of neighbours to compare (kept within 1 and n-1).
#'
#' @return A numeric vector of agreement scores between 0 and 1; nodes
#'   with NA coordinates score NA.#'
#' @keywords internal
.knnAgreement <- function(codes, xy, k) {
    codes <- as.matrix(codes)
    xy <- as.matrix(xy)
    k <- max(1L, min(as.integer(k), nrow(codes) - 1L))
    fullNN <- .knnIndices(codes, k)
    xyNN <- .knnIndices(xy, k)
    vapply(seq_len(nrow(codes)), function(i) {
        if (anyNA(xy[i, ])) {
            return(NA_real_)
        }
        a <- fullNN[[i]]
        b <- xyNN[[i]]
        if (length(a) == 0L && length(b) == 0L) {
            return(1) # degenerate: no neighbours to compare in either space
        }
        length(intersect(a, b)) / length(union(a, b))
    }, numeric(1))
}

#' k Nearest Neighbour Indices of Each Row
#'
#' @param m Numeric matrix.
#' @param k Number of neighbours to return per row.
#'
#' @return A list with one integer vector of neighbour row indices per row.
#'
#' @keywords internal
.knnIndices <- function(m, k) {
    d <- stats::dist(m)
    if (nrow(m) < 2L) {
        return(rep(list(integer(0)), nrow(m)))
    }
    # order() is stable for equal keys, so ties break by row index. Ties
    # that differ only in floating-point noise across the two spaces are
    # quantized away by comparing rounded distances.
    ord <- apply(as.matrix(round(d, 10L)), 1L, order, method = "radix")
    lapply(seq_len(nrow(m)), function(i) {
        others <- ord[, i]
        others <- others[others != i]
        others[seq_len(min(k, length(others)))]
    })
}

drawProjection <- function(df, rs, colorbyGroups, sce, outputList = list(),
    agreement = NULL) {
    colN <- names(df)[seq_len(2)]
    df <- .projectionNodeData(df, sce)
    grp <- .projectionGroups(df, colorbyGroups, outputList)
    df <- grp$df
    if (!is.null(agreement)) {
        # agreement is indexed by node id; .projectionNodeData ordered by id.
        df$agreement <- agreement[df$id]
    }
    aesArgs <- list(
        x = rlang::quo(.data[[colN[1]]]), y = rlang::quo(.data[[colN[2]]]),
        # key identifies the node for plotly selections.
        key = rlang::quo(.data$id),
        text = rlang::quo({
            agreementLine <- if (!is.null(agreement)) {
                paste("<br>agreement:", format(.data$agreement, digits = 2))
            } else {
                ""
            }
            paste(
                "cluster:", .data$cluster, "<br>N cells:", .data$N,
                "<br>mean:", format(.data$mean, digits = 2),
                "<br>3rd Q:", format(.data$thrdQu, digits = 2),
                "<br>max:", format(.data$max, digits = 2),
                agreementLine
            )
        })
    )
    if (!is.null(agreement) && grp$nGrps == 1L) {
        aesArgs$color <- rlang::quo(.data$agreement)
        p3 <- ggplot2::ggplot(df, do.call(ggplot2::aes, aesArgs)) +
            ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
            ggplot2::scale_colour_viridis_c(name = "agreement", direction = -1)
    } else if (grp$nGrps == 1L) {
        p3 <- ggplot2::ggplot(df, do.call(ggplot2::aes, aesArgs)) +
            ggplot2::geom_point(color = "lightblue", size = 1.5,
                show.legend = FALSE, na.rm = TRUE)
    } else {
        aesArgs$color <- rlang::quo(.data$colGrp)
        p3 <- ggplot2::ggplot(df, do.call(ggplot2::aes, aesArgs)) +
            ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
            ggplot2::scale_color_manual(values = grp$colors, drop = FALSE)
    }
    # Red overlay of the selected nodes, purely visual.
    selected <- df[df$id %in% as.integer(rs), , drop = FALSE]
    if (nrow(selected) > 0L) {
        p3 <- p3 + ggplot2::geom_point(
            data = selected,
            mapping = ggplot2::aes(x = .data[[colN[1]]], y = .data[[colN[2]]]),
            color = "red", size = 0.8, inherit.aes = FALSE, na.rm = TRUE
        )
    }
    p3
}

#' Per-Node Data of a Projection Plot
#'
#' Joins the SOM node statistics when \code{df} has no \code{id}, orders the
#' rows by node (row i is node i, as the plotly point-number fallback
#' expects) and adds the column aliases used in the hover text.
#'
#' @param df Projection coordinates (first two columns), optionally with
#'   \code{id} and node statistics.
#' @param sce \code{SingleCellExperiment} with \code{SOM_stats} metadata.
#'
#' @return The completed data frame with \code{cluster}, \code{N} and
#'   \code{thrdQu}.
#'
#' @keywords internal
.projectionNodeData <- function(df, sce) {
    if (!"id" %in% names(df)) {
        df$id <- seq_len(nrow(df))
        df <- dplyr::left_join(df, S4Vectors::metadata(sce)$SOM_stats,
            by = "id")
    }
    df <- df[order(df$id), ]
    if ("n" %in% names(df) && !"N" %in% names(df)) df$N <- df$n
    if ("rdQu" %in% names(df) && !"thrdQu" %in% names(df)) {
        df$thrdQu <- df$rdQu
    }
    if ("thirdQu" %in% names(df) && !"thrdQu" %in% names(df)) {
        df$thrdQu <- df$thirdQu
    }
    df$cluster <- df$id
    df
}

#' Colour Groups of a Projection Plot
#'
#' @param df Per-node data with \code{cluster}.
#' @param colorbyGroups Saved group names to colour by.
#' @param outputList Named list of saved groups.
#'
#' @return A list with \code{df} (with factor \code{colGrp} when grouping;
#'   nodes in several groups get a combined label), the number of groups
#'   \code{nGrps} and the palette \code{colors}.
#'
#' @keywords internal
.projectionGroups <- function(df, colorbyGroups, outputList) {
    if (length(colorbyGroups) < 1L) {
        return(list(df = df, nGrps = 1L, colors = NULL))
    }
    df$colGrp <- "other"
    for (cg in colorbyGroups) {
        mask <- df$cluster %in% outputList[[cg]]
        df$colGrp[mask] <- ifelse(df$colGrp[mask] == "other", cg,
            paste(df$colGrp[mask], cg))
    }
    df$colGrp <- factor(df$colGrp)
    nGrps <- length(levels(df$colGrp))
    list(df = df, nGrps = nGrps,
        colors = colorRampPalette(brewer.pal(8, "Set2"))(nGrps + 1L))
}


plotViolin2Func <- function(sce, somCodesName = "SOM_codes", violinSelection,
    upsetSelection, outputList) {
    if (length(upsetSelection) < 3) {
        upsetSelection <- names(outputList)
    }
    data <- .violinData(sce@metadata[[somCodesName]], upsetSelection,
        outputList, violinSelection)
    .violinGgplot(data, "No groups with \u2265 2 SOM nodes selected",
        xVar = "grpName")
}


upsetPlotFunc <- function(upsetSelection, outputList, sce, maxCombs = 100) {
    if (length(upsetSelection) < 2) {
        return(NULL)
    }
    if (is.null(grDevices::dev.list())) {
        grDevices::pdf(NULL)
        on.exit(grDevices::dev.off(), add = TRUE)
    }
    if (length(upsetSelection) > 31) upsetSelection <- head(upsetSelection, 31)
    cm <- ComplexHeatmap::make_comb_mat(outputList[upsetSelection])
    nCombs <- length(ComplexHeatmap::comb_name(cm))
    if (nCombs > maxCombs) {
        stop(
            "Too many distinct group overlap combinations (", nCombs,
            ") to render an UpSet plot; select fewer/less-overlapping groups ",
            "(currently allows up to ", maxCombs, ")."
        )
    }
    ann <- .upsetAnnotations(cm, upsetSelection, outputList,
        S4Vectors::metadata(sce)$SOM_stats)
    ComplexHeatmap::UpSet(cm,
        comb_order = order(
            ComplexHeatmap::comb_degree(cm),
            -ComplexHeatmap::comb_size(cm)
        ),
        top_annotation = ann$top,
        right_annotation = ann$side
    )
}

#' Cell and Node Count Annotations of the UpSet Plot
#'
#' @param cm Combination matrix (\code{ComplexHeatmap::make_comb_mat()}).
#' @param upsetSelection Group names shown.
#' @param outputList Named list of SOM node ids per group.
#' @param somStats SOM node statistics with cell counts \code{n}.
#'
#' @return A list with the \code{top} (per combination) and \code{side}
#'   (per group) annotations.
#'
#' @keywords internal
.upsetAnnotations <- function(cm, upsetSelection, outputList, somStats) {
    cellsOf <- function(nodes) sum(somStats[as.integer(nodes), "n"])
    combs <- ComplexHeatmap::comb_name(cm)
    ncells <- vapply(combs, function(cn) {
        cellsOf(ComplexHeatmap::extract_comb(cm, cn))
    }, numeric(1))
    grpCells <- vapply(upsetSelection, function(g) cellsOf(outputList[[g]]),
        numeric(1))
    grpSoms <- unlist(lapply(outputList[upsetSelection], length))
    bar <- function(v, fill, width) {
        ComplexHeatmap::anno_barplot(v, add_numbers = TRUE,
            gp = grid::gpar(fill = fill), width = grid::unit(width, "cm"))
    }
    list(
        top = ComplexHeatmap::HeatmapAnnotation(
            "cell #" = bar(ncells, "black", 4),
            "som #" = bar(ComplexHeatmap::comb_size(cm), "black", 4),
            gap = grid::unit(2, "mm"), annotation_name_side = "left",
            annotation_name_rot = 0
        ),
        side = ComplexHeatmap::rowAnnotation(
            "som #" = bar(grpSoms, NULL, 3),
            "cell #" = bar(grpCells, NULL, 3),
            gap = grid::unit(2, "mm"), annotation_name_rot = 0
        )
    )
}
