#' @rdname plotSOMScatter
#' @title Scatter plot
#'
#' @description Bivariate scatter plots including visualization of
#' (group-specific) gates, their boundaries and percentage of selected cells.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character vector specifying which channels to plot.
#' @param metaSlot name of the metadata slot containing SOM codes.
#' @param pointSize column in SOM stats used to size points.
#' @param color_by column to color points by.
#' @param bins number of bins when coloring by density.
#' @param assay name of the assay to use.
#' @param statsSlot name of the metadata slot containing SOM stats.
#' @param label how axis labels should be constructed.
#' @param zeros logical specifying whether to include 0 values.
#' @param k_pal optional cluster color palette.
#' @param xRN optional row names.
#' @param xCN optional channel names.
#'
#' @return a \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#' @importFrom ggplot2 scale_color_manual facet_wrap facet_grid vars
#' @importFrom ggplot2 guides guide_legend ylab theme_bw theme
#' @importFrom ggplot2 element_blank element_text element_rect
#' @importFrom grid unit
#' @importFrom rlang sym .data
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom SingleCellExperiment int_colData
#' @importFrom CATALYST channels cluster_codes cluster_ids
#' @importFrom S4Vectors metadata
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' sce <- CySA_example_sce()
#' plotSOMScatter(sce, chs = c("marker1", "marker2"))
#'
#' @export
plotSOMScatter <- function(x, chs, metaSlot = "SOM_codes", pointSize = "n",
    color_by = "n", bins = 100, assay = "exprs", statsSlot = "SOM_stats",
    label = c("target", "channel", "both"), zeros = FALSE, k_pal = NULL,
    xRN = NULL, xCN = NULL) {
    label <- match.arg(label)
    pointSize <- pointSize[1]
    color_by <- color_by[1]
    if (is.null(xRN)) xRN <- rownames(x)
    if (is.null(xCN)) xCN <- CATALYST::channels(x)
    inp <- .somScatterInputs(x, chs, metaSlot, statsSlot, xRN, xCN, assay)
    if (is.null(inp)) {
        return(NULL)
    }
    dat <- .somScatterData(inp, label, xRN, xCN, zeros)
    col <- .somScatterColour(dat$df, color_by, x, k_pal)
    if (!pointSize %in% colnames(dat$df)) pointSize <- NULL
    aesArgs <- list(x = rlang::sym(dat$chs[1]), y = rlang::sym(dat$chs[2]))
    if (!is.null(pointSize)) aesArgs$size <- rlang::sym(pointSize)
    if (!is.null(col$var)) aesArgs$color <- rlang::sym(col$var)
    p1 <- ggplot2::ggplot(dat$df, do.call(ggplot2::aes, aesArgs)) +
        ggplot2::geom_point(alpha = 0.2, na.rm = TRUE)
    for (layer in list(col$scale, col$guides, dat$ylab, dat$facet)) {
        if (!is.null(layer)) p1 <- p1 + layer
    }
    p1 + ggplot2::theme_bw() + ggplot2::theme(
        aspect.ratio = 1,
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(color = "black"),
        strip.background = ggplot2::element_rect(fill = "white"),
        legend.key.height = grid::unit(0.8, "lines")
    )
}

#' Validated Inputs of plotSOMScatter()
#'
#' Checks the metadata slots and channels; channels missing from the SOM
#' codes but present in the assay get per-node mean codes computed locally
#' (\code{x} is not modified).
#'
#' @inheritParams plotSOMScatter
#'
#' @return \code{NULL} (with a warning) when no requested channel is known,
#'   otherwise a list with \code{codes} (SOM codes incl. computed columns),
#'   \code{stats} (SOM stats or \code{NULL}) and \code{chs} (known channels).
#'
#' @keywords internal
.somScatterInputs <- function(x, chs, metaSlot, statsSlot, xRN, xCN, assay) {
    md <- S4Vectors::metadata(x)
    if (!metaSlot %in% names(md)) {
        stop("Need '", metaSlot, "' in metadata of x", call. = FALSE)
    }
    codes <- md[[metaSlot]]
    if (!is.null(statsSlot) && !statsSlot %in% names(md)) {
        warning("'", statsSlot,
            "' not found in metadata of x; proceeding without stats")
        statsSlot <- NULL
    }
    missing <- setdiff(chs, unique(c(colnames(codes), xRN, xCN)))
    if (length(missing) == length(chs)) {
        warning(
            "Unknown channels/markers: ", paste(missing, collapse = ", "),
            "\nAvailable in SOM_codes: ",
            paste(colnames(codes), collapse = ", "),
            call. = FALSE
        )
        return(NULL)
    }
    if (length(missing)) {
        warning("Unknown channels/markers dropped: ",
            paste(missing, collapse = ", "), call. = FALSE)
        chs <- setdiff(chs, missing)
    }
    y <- SummarizedExperiment::assay(x, assay)
    cid <- SummarizedExperiment::colData(x)$cluster_id
    for (ch in setdiff(chs, colnames(codes))) {
        if (!ch %in% rownames(y)) {
            stop("Cannot compute SOM stats for '", ch,
                "': not found in assay rows.", call. = FALSE)
        }
        warning("computing SOM stats for '", ch, "'")
        vals <- vapply(as.integer(unique(cid)), function(cl) {
            mean(y[ch, which(cid == cl), drop = TRUE])
        }, numeric(1))
        codes <- cbind(codes, stats::setNames(vals, NULL))
        colnames(codes)[ncol(codes)] <- ch
    }
    list(codes = codes, chs = chs, stats = .somScatterStats(x, md, statsSlot))
}

#' Stats Frame of plotSOMScatter()
#'
#' The SOM statistics joined with the known-label node summary when the
#' \code{node_labels} slot exists, so \code{purity} and
#' \code{cell_type_majority} can drive point colour and size.
#'
#' @param x The \code{SingleCellExperiment} being plotted.
#' @param md Its metadata (as a list).
#' @param statsSlot Name of the SOM statistics metadata slot.
#'
#' @return The joined statistics data frame, or \code{NULL} when no
#'   statistics are available.
#'
#' @keywords internal
.somScatterStats <- function(x, md, statsSlot) {
    stats <- if (!is.null(statsSlot)) md[[statsSlot]]
    if (is.null(stats) || is.null(md$node_labels)) {
        return(stats)
    }
    dplyr::left_join(stats, md$node_labels, by = "id")
}

#' Plot Data of plotSOMScatter()
#'
#' Relabels the SOM code columns, joins the node statistics, drops nodes
#' with zeros (unless \code{zeros}) and, for more than two channels, melts
#' the data for faceting.
#'
#' @param inp Output of \code{.somScatterInputs()}.
#' @inheritParams plotSOMScatter
#'
#' @return A list with \code{df}, the plotted channel names \code{chs} and
#'   the optional \code{facet} and \code{ylab} layers.
#'
#' @keywords internal
.somScatterData <- function(inp, label, xRN, xCN, zeros) {
    chs <- inp$chs
    yy <- inp$codes[, chs, drop = FALSE]
    nms <- switch(label,
        target = xRN,
        channel = xCN,
        both = ifelse(xCN == xRN, xCN, paste(xCN, xRN, sep = "-"))
    )
    # Rename only the columns that correspond to rownames(x) / channels.
    renameMap <- stats::setNames(nms, xRN)
    found <- chs %in% names(renameMap)
    colnames(yy)[found] <- renameMap[chs[found]]
    if (!is.null(inp$stats)) {
        df <- cbind(yy, inp$stats)
    } else {
        df <- as.data.frame(yy)
        df$n <- 1
    }
    if (is.null(chs) || length(chs) < 2) {
        stop("At least two channels are required for a scatter plot.",
            call. = FALSE)
    }
    dropZeros <- function(df, chs) {
        if (zeros) {
            return(df)
        }
        df[rowSums(df[, chs[c(1, 2)], drop = FALSE] == 0) == 0, ,
            drop = FALSE]
    }
    df <- dropZeros(df, chs)
    facet <- ylab <- NULL
    if (length(chs) > 2) {
        idVars <- unique(c(chs[1], colnames(inp$stats)))
        df <- reshape2::melt(df, id.vars = idVars)
        facet <- ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym("variable")))
        ylab <- ggplot2::ylab(NULL)
        chs[2] <- "value"
        df <- dropZeros(df, chs)
    }
    list(df = df, chs = chs, facet = facet, ylab = ylab)
}

#' Colour Scale of plotSOMScatter()
#'
#' @param df Plot data.
#' @inheritParams plotSOMScatter
#'
#' @return A list with the colour column \code{var} (or \code{NULL}), the
#'   colour \code{scale} and the legend \code{guides} (or \code{NULL}).
#'
#' @keywords internal
.somScatterColour <- function(df, color_by, x, k_pal) {
    colVar <- if (color_by %in% colnames(df)) color_by else NULL
    out <- list(var = colVar, scale = NULL, guides = NULL)
    if (is.null(colVar)) {
        return(out)
    }
    if (is.numeric(df[[colVar]])) {
        out$scale <- ggplot2::scale_color_gradientn(
            colors = c("navy", rev(RColorBrewer::brewer.pal(11, "Spectral"))),
            guide = "colorbar"
        )
    } else if (colVar %in% names(CATALYST::cluster_codes(x))) {
        if (is.null(k_pal)) k_pal <- CySA_default_cluster_cols()
        out$scale <- ggplot2::scale_color_manual(values = k_pal)
        out$guides <- ggplot2::guides(col = ggplot2::guide_legend(
            override.aes = list(alpha = 1, size = 3)
        ))
    } else {
        # Categorical property column (e.g. the known-label majority): map
        # every level to a Set2 colour, cycled when there are more levels
        # than the palette holds; unknown values render grey.
        levs <- unique(stats::na.omit(df[[colVar]]))
        pal <- RColorBrewer::brewer.pal(max(3, min(8, length(levs))), "Set2")
        out$scale <- ggplot2::scale_color_manual(
            values = stats::setNames(
                rep(pal, length.out = length(levs)), levs),
            na.value = "grey70"
        )
        out$guides <- ggplot2::guides(col = ggplot2::guide_legend(
            override.aes = list(alpha = 1, size = 3)
        ))
    }
    out
}
