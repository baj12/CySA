# Node-properties 2D plot: one point per SOM node, axes from the per-node
# property table (SOM_stats + optional known-label summary + ambiguity).

#' Node-Properties 2D Plot of the SOM Nodes
#'
#' Scatter plot in which each point is a SOM node, plotted against two of
#' its per-node properties. The default axes pair the known-label purity
#' with the node size when labels are available, and the per-node ambiguity
#' (median own/rival distance ratio) with the node size otherwise, so the
#' plot works with or without annotations. Points carry the node id as
#' their plotly \code{key}: clicks and lasso selections made in the Shiny
#' app resolve to nodes through the same selection pipeline as the other 2D
#' plots.
#'
#' @param x A \code{SingleCellExperiment} whose metadata holds the SOM
#'   codes and \code{SOM_stats}; \code{node_labels} (added by
#'   \code{\link{prepClusterSelectorData}} when \code{cellTypeCol} is
#'   given) is used when present.
#' @param xVar,yVar Names of numeric property columns for the axes;
#'   \code{NULL} resolves to \code{"purity"} x \code{"n"} when \code{purity}
#'   is available and \code{"medianRatio"} x \code{"n"} otherwise.
#' @param somCodesName Name of the SOM codes metadata slot (used by the
#'   on-the-spot ambiguity computation).
#' @param somStatsSlot Name of the per-node statistics metadata slot.
#' @param labelSlot Name of the known-label node summary metadata slot;
#'   joined when present.
#' @param ambiguity Optional precomputed per-node ambiguity (the output of
#'   \code{\link{nodeAmbiguity}}). When \code{NULL} it is computed from
#'   \code{x} on the spot, which is potentially expensive on large objects;
#'   pass a subsample's summary or a precomputed one in that case.
#' @param ambiguousThreshold Ambiguity threshold, as in
#'   \code{\link{nodeAmbiguity}}.
#' @param maxCells Cell cap for the on-the-spot ambiguity computation.
#' @param selectedIds Integer node ids to highlight in red (the current
#'   selection in the Shiny app).
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{nodeAmbiguity}},
#'   \code{\link{prepClusterSelectorData}}
#'
#' @export
plotNodeProps <- function(x, xVar = NULL, yVar = NULL,
    somCodesName = "SOM_codes", somStatsSlot = "SOM_stats",
    labelSlot = "node_labels", ambiguity = NULL,
    ambiguousThreshold = 0.9, maxCells = 20000, selectedIds = NULL) {
    stats <- S4Vectors::metadata(x)[[somStatsSlot]]
    if (is.null(stats)) {
        stop("No '", somStatsSlot, "' found in metadata(x)")
    }
    if (!"id" %in% colnames(stats)) {
        stop("'", somStatsSlot, "' has no 'id' column")
    }
    if (is.null(ambiguity)) {
        ambiguity <- nodeAmbiguity(x, somCodesName = somCodesName,
            ambiguousThreshold = ambiguousThreshold, maxCells = maxCells)
    }
    labels <- S4Vectors::metadata(x)[[labelSlot]]
    df <- stats
    if (!is.null(labels)) {
        df <- left_join(df, labels, by = "id")
    }
    # ambiguity also has an 'n' (cells sampled); the plot's 'n' is the node
    # size from SOM_stats, so keep the left-hand name on collision.
    df <- left_join(df, ambiguity, by = "id", suffix = c("", ".amb"))
    numCols <- names(df)[vapply(df, is.numeric, logical(1))]
    if (is.null(xVar)) {
        xVar <- if ("purity" %in% numCols && any(!is.na(df$purity))) {
            "purity"
        } else {
            "medianRatio"
        }
    }
    if (is.null(yVar)) {
        yVar <- "n"
    }
    unknown <- setdiff(c(xVar, yVar), numCols)
    if (length(unknown) > 0L) {
        stop("Unknown or non-numeric property column(s): ",
            paste(unknown, collapse = ", "))
    }
    hover <- paste0("node: ", df$id, "<br>n cells: ", df$n)
    if ("purity" %in% colnames(df)) {
        hover <- paste0(hover,
            "<br>purity: ", format(df$purity, digits = 2))
    }
    if ("cell_type_majority" %in% colnames(df)) {
        hover <- paste0(hover,
            "<br>majority: ", df$cell_type_majority)
    }
    df$text <- hover
    p <- ggplot(df, aes(x = .data[[xVar]], y = .data[[yVar]],
        key = .data$id, text = .data$text))
    if ("purity" %in% colnames(df)) {
        p <- p +
            geom_point(mapping = aes(colour = .data$purity), size = 1.5,
                na.rm = TRUE) +
            scale_colour_viridis_c(name = "purity", direction = -1,
                na.value = "grey70")
    } else {
        p <- p +
            geom_point(color = "lightblue", size = 1.5,
                show.legend = FALSE, na.rm = TRUE)
    }
    selected <- df[df$id %in% as.integer(selectedIds), , drop = FALSE]
    if (nrow(selected) > 0L) {
        p <- p +
            geom_point(data = selected,
                mapping = aes(x = .data[[xVar]], y = .data[[yVar]]),
                color = "red", size = 0.8, inherit.aes = FALSE,
                na.rm = TRUE)
    }
    p + labs(x = xVar, y = yVar) +
        theme_bw() +
        theme(
            aspect.ratio = 1,
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white"),
            legend.key.height = grid::unit(0.8, "lines")
        )
}
