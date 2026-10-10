# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# geom-ridgeline-traced.R ----
# A ggridges-derived ridgeline geom that additionally draws connected
# "trace" polylines and outlier dots inside the same layer: ggridges'
# GeomRidgeline splits one layer's data by the "datatype" column
# (ridgeline / vline / point rows, each rendered by its own grob maker);
# here a fourth datatype, "path", is understood and rendered as polylines
# across the ridges. Rows of the path datatype carry x, y (the vertex
# position, already offset above the ridge base) and group (one group per
# polyline); their ymin is set below every ridge base so the trace group
# sorts last in GeomRidgeline's draw_panel ordering and is drawn on top.

# ggname: grob naming helper, as in ggridges (internal there, so a local
# copy from ggplot2's utils is kept here).
ggname <- function(prefix, grob) {
    grob$name <- grid::grobName(grob, prefix)
    grob
}

#' Ridgeline Geom with Connected Trace Polylines
#'
#' \code{geom_ridgeline_traced()} behaves like
#' \code{ggridges::geom_ridgeline()} but understands one additional value of
#' the \code{datatype} aesthetic: rows with \code{datatype = "path"} are
#' drawn as polylines connecting their vertices (one polyline per
#' \code{group}), drawn above the ridge lines; rows with
#' \code{datatype = "point"} are drawn as the outlier dots of the traces,
#' like in \code{ggridges}. \code{GeomRidgelineTraced} is the underlying
#' ggproto object, a subclass of \code{ggridges::GeomRidgeline}.
#'
#' Rows of the \code{"path"} datatype should carry \code{x}, \code{y} (the
#' vertex position, already offset above the ridge base) and \code{group}
#' (one group per polyline). Rows with the \code{"point"} datatype need
#' \code{x} and \code{y} (the dot position).
#'
#' @param mapping Aesthetic mappings created by \code{ggplot2::aes()}.
#' @param data Layer data; \code{NULL} inherits the plot data.
#' @param stat,position Statistical transformation and position
#'   adjustment (both \code{"identity"}).
#' @param na.rm Logical; drop rows with missing values silently.
#' @param show.legend Logical or \code{NA}; include the layer in legends.
#' @param inherit.aes Logical; combine with the plot's default aesthetics.
#' @param ... Further parameters passed to \code{ggplot2::layer()}.
#'
#' @return A ggplot2 layer.
#'
#' @keywords internal
#' @rdname geom_ridgeline_traced
geom_ridgeline_traced <- function(mapping = NULL, data = NULL,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE, ...) {
    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = ggplot2::StatIdentity,
        geom = GeomRidgelineTraced,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}

#' @rdname geom_ridgeline_traced
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 ggproto Geom layer StatIdentity alpha coord_munch
#'     fill_alpha
#' @importFrom grid gpar grobName grobTree nullGrob pointsGrob polylineGrob
#' @importFrom ggridges GeomRidgeline
GeomRidgelineTraced <- ggplot2::ggproto("GeomRidgelineTraced",
    ggridges::GeomRidgeline,
    # Metadata columns (marker/cluster/cellId/isExtreme) ride along as
    # aesthetics so they survive into the build data, where tests and
    # grob makers read them; NULL defaults keep them data-driven only.
    # The parent default_aes is merged in first (ggproto concatenation).
    default_aes = ggplot2::aes(
        marker = NULL, cluster = NULL, cellId = NULL, isExtreme = NULL,
        y = 0, colour = "black", fill = "grey70", linewidth = 0.5,
        linetype = 1, min_height = 0, scale = 1, alpha = NA,
        datatype = "ridgeline",
        point_shape = 19, point_size = 1.5, point_stroke = 0.5,
        point_colour = NULL, point_fill = NULL, point_alpha = NULL,
        vline_colour = NULL, vline_width = NULL, vline_linetype = NULL,
        vline_size = NULL
    ),
    setup_data = function(self, data, params) {
        if (!"scale" %in% names(data)) {
            data$scale <- 1
        }
        if (!"min_height" %in% names(data)) {
            data$min_height <- 0
        }
        if (!"ymin" %in% names(data)) {
            data$ymin <- NA_real_
        }
        if (!"ymax" %in% names(data)) {
            data$ymax <- NA_real_
        }
        isTrace <- !is.na(data$datatype) &
            data$datatype %in% c("path", "point")
        # Ridge rows take the ggridges identity encoding (y is the ridge
        # base, height stacks upward); trace rows keep the ymin they carry
        # for draw ordering (below every ridge base, so the draw_panel
        # ordering draws them last, on top) and their absolute vertex
        # position in y.
        data$ymin[!isTrace] <- data$y[!isTrace]
        data$ymax[!isTrace] <- data$y[!isTrace] +
            data$scale[!isTrace] * data$height[!isTrace]
        if (any(isTrace)) {
            if (any(is.na(data$ymin[isTrace]))) {
                # ymin is a position aesthetic and trains the y scale, so
                # the draw-order key must sit only a hair below the lowest
                # ridge base: a fixed offset (formerly 1 data unit) dwarfs
                # small density values and squashes every ridge into a
                # flat line at the top of the panel.
                yAll <- c(data$ymin[!isTrace], data$y)
                lo <- min(0, yAll, na.rm = TRUE)
                span <- diff(range(c(lo, yAll, data$ymax[!isTrace]),
                    na.rm = TRUE))
                if (!is.finite(span) || span <= 0) {
                    span <- 1
                }
                data$ymin[isTrace & is.na(data$ymin)] <- lo - 1e-3 * span
            }
            if (any(is.na(data$ymax[isTrace]))) {
                data$ymax[isTrace & is.na(data$ymax)] <-
                    data$y[isTrace & is.na(data$ymax)]
            }
        }
        data
    },

    draw_group = function(self, data, panel_params, coord, na.rm = FALSE) {
        if (na.rm) {
            data <- data[stats::complete.cases(data[c("x", "ymin", "ymax")]), ]
        }
        data_list <- split(data, factor(data$datatype))
        point_grob <- self$make_point_grob(data_list[["point"]],
            panel_params, coord)
        vline_grob <- self$make_vline_grob(data_list[["vline"]],
            panel_params, coord)
        path_grob <- self$make_path_grob(data_list[["path"]], panel_params,
            coord)
        ridges <- data_list[["ridgeline"]]
        if (is.null(ridges)) {
            return(grid::grobTree(vline_grob, path_grob, point_grob))
        }
        rg <- self$make_ridge_grobs(ridges, panel_params, coord)
        # Traces above the ridge lines, dots above the traces.
        grid::grobTree(rg$area, vline_grob, rg$line, path_grob, point_grob)
    },

    # Connected traces: one polyline per group, vertices ordered by x.
    make_path_grob = function(data, panel_params, coord) {
        if (is.null(data) || nrow(data) == 0L) {
            return(grid::nullGrob())
        }
        data <- data[order(data$group, data$x), , drop = FALSE]
        path_aes <- unique(data[c("colour", "linewidth", "linetype", "alpha")])
        if (nrow(path_aes) > 1L) {
            stop("Aesthetics can not vary along a trace polyline")
        }
        munched <- ggplot2::coord_munch(coord, data, panel_params)
        ggname("geom_ridgeline_traced", grid::polylineGrob(
            munched$x, munched$y, id = munched$group,
            default.units = "native",
            gp = grid::gpar(
                col = ggplot2::alpha(path_aes$colour, path_aes$alpha),
                lwd = path_aes$linewidth * .pt,
                lty = path_aes$linetype
            )
        ))
    },

    # Area and outline of one ridge (block adapted from
    # ggridges::GeomRidgeline$draw_group).
    make_ridge_grobs = function(self, data, panel_params, coord) {
        data <- data[order(data$group, data$x), ]
        # Remove all points that fall below the minimum height.
        data$ymax[data$height < data$min_height] <- NA
        aes <- unique(data[c("colour", "fill", "linewidth", "linetype",
            "alpha")])
        if (nrow(aes) > 1) {
            stop("Aesthetics can not vary along a ridgeline")
        }
        # Instead of removing NA values and plotting a single polygon, stop
        # the polygon at missing values and start a new one after them.
        missing_pos <- !stats::complete.cases(data[c("x", "ymin", "ymax")])
        ids <- cumsum(missing_pos) + 1
        ids[missing_pos] <- NA
        munched_poly <- ggplot2::coord_munch(coord, data.frame(
            x = c(data$x, rev(data$x)),
            y = c(data$ymax, rev(data$ymin)),
            id = c(ids, rev(ids))
        ), panel_params)
        munched_line <- ggplot2::coord_munch(coord,
            data.frame(x = data$x, y = data$ymax, id = ids), panel_params)
        list(
            area = self$make_area_grob(munched_poly, as.list(aes)),
            line = self$make_line_grob(munched_line, munched_poly,
                as.list(aes))
        )
    },

    make_point_grob = function(data, panel_params, coord) {
        if (is.null(data)) {
            return(grid::nullGrob())
        }
        # Trace dots sit at their vertex position (y), not at ymin: the
        # ymin column of trace rows carries the draw-order key, so the
        # ggridges default (data$y <- data$ymin) would park the dots
        # below the panel.
        coords <- coord$transform(data, panel_params)
        pointColour <- if (is.null(data$point_colour)) {
            if (is.null(data$point_color)) {
                data$colour
            } else {
                data$point_color
            }
        } else {
            data$point_colour
        }
        pointFill <- if (is.null(data$point_fill)) {
            data$fill
        } else {
            data$point_fill
        }
        pointAlpha <- if (is.null(data$point_alpha)) {
            data$alpha
        } else {
            data$point_alpha
        }
        ggname("geom_ridgeline_traced", grid::pointsGrob(
            coords$x, coords$y,
            pch = coords$point_shape,
            gp = grid::gpar(
                col = ggplot2::alpha(pointColour, pointAlpha),
                fill = ggplot2::fill_alpha(pointFill, pointAlpha),
                fontsize = coords$point_size * .pt +
                    coords$point_stroke * .stroke / 2,
                lwd = coords$point_stroke * .stroke / 2
            )
        ))
    }
)
