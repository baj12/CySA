# tests/testthat/test-geom-ridgeline-traced.R
#
# Tests for the ggridges-derived ridgeline geom that additionally draws
# connected outlier traces: one data frame with mixed "datatype" rows
# (ridgeline / path / point) renders ridges, polylines across ridges and
# outlier dots in a single layer, reusing ggridges' draw machinery.

library(testthat)

# ---- Tiny hand-checkable fixture --------------------------------------------

# Three ridges on a shared x range, plus a two-vertex trace whose first
# vertex sits on ridge A (base 0 + density 1 = y 1) and second on ridge B
# (base 2 + density 1 = y 3); the plot function pre-computes vertex
# positions, the geom draws them where the data puts them.
make_layer_frame <- function() {
    data.frame(
        x = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 2),
        y = c(0, 0, 0, 2, 2, 2, 4, 4, 4, 1, 3),
        # density heights on the shared grid; the trace sits at x = 0 (ridge
        # A: height 1) and x = 2 (ridge B: height 1).
        height = c(1, 0.5, 1, 1, 0.5, 1, 1, 0.5, 1, NA, NA),
        datatype = c(rep("ridgeline", 9), rep("path", 2)),
        group = c(rep(1L, 3), rep(2L, 3), rep(3L, 3), rep(10L, 2)),
        colour = "red",
        fill = "grey80",
        stringsAsFactors = FALSE
    )
}

test_that("GeomRidgelineTraced inherits from ggridges GeomRidgeline", {
    expect_s3_class(ggridges::GeomRidgeline, "ggproto")
    expect_identical(
        CySA:::GeomRidgelineTraced$super(),
        ggridges::GeomRidgeline
    )
})

test_that("path rows draw a polyline across ridges in one layer", {
    d <- make_layer_frame()
    p <- ggplot2::ggplot(d) +
        ggplot2::aes(x = x, y = y, height = height, datatype = datatype) +
        geom_ridgeline_traced()

    b <- ggplot2::ggplot_build(p)
    # Exactly one layer carries all three datatypes: the ggridges stat is
    # not needed here, so a single GeomRidgelineTraced layer holds ridges,
    # paths and points.
    expect_length(p$layers, 1)
    expect_identical(class(p$layers[[1]]$geom)[1], "GeomRidgelineTraced")

    ld <- b$data[[1]]
    expect_setequal(unique(ld$datatype), c("ridgeline", "path"))
    # Path group must span the ridge rows: two vertices with different y.
    pathRows <- ld[ld$datatype == "path", ]
    expect_equal(nrow(pathRows), 2)
    expect_false(isTRUE(all.equal(pathRows$y[1], pathRows$y[2])))
})

test_that("path vertices land on ridge tops at interpolated density", {
    d <- make_layer_frame()
    p <- ggplot2::ggplot(d) +
        ggplot2::aes(x = x, y = y, height = height, datatype = datatype) +
        geom_ridgeline_traced()

    b <- ggplot2::ggplot_build(p)
    # Ridge 1 spans x = 0..2 with height 1 at x = 0; a path vertex at x = 0
    # must sit at y = base + density = 0 + 1.
    ld <- b$data[[1]]
    v1 <- ld[ld$datatype == "path" & ld$x == 0, ]
    expect_equal(v1$y, 1)
    # Ridge 2 spans x = 0..2 with height 1 at x = 2 (base 2): vertex at
    # y = 2 + 1 = 3.
    v2 <- ld[ld$datatype == "path" & ld$x == 2, ]
    expect_equal(v2$y, 3)
})

test_that("point rows draw dots with default point aesthetics", {
    d <- make_layer_frame()
    d$datatype <- c(rep("ridgeline", 9), "path", "point")
    p <- ggplot2::ggplot(d) +
        ggplot2::aes(x = x, y = y, height = height, datatype = datatype) +
        geom_ridgeline_traced()

    b <- ggplot2::ggplot_build(p)
    expect_identical(class(p$layers[[1]]$geom)[1], "GeomRidgelineTraced")
    # The build data keeps the point row.
    ld <- b$data[[1]]
    expect_true("point" %in% ld$datatype)
    expect_equal(nrow(ld[ld$datatype == "point", ]), 1)
})

test_that("geom_ridgeline_traced is exported and documented", {
    # NAMESPACE must export the new geom and ggproto so users and tests can
    # call it; roxygen2 generates both from the R source.
    ns <- asNamespace("CySA")
    expect_true(exists("geom_ridgeline_traced", ns, inherits = FALSE))
    expect_true(exists("GeomRidgelineTraced", ns, inherits = FALSE))
})
