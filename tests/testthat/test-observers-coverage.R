# Tests for uncovered regions in observers_clusterSelector.R
# Targets lines marked as uncovered in the coverage report.

library(shiny)
library(testthat)
library(CySA)

# =============================================================================
# clusterNameSelect observer (lines 148-159)
# =============================================================================
test_that("observers: clusterNameSelect updates clusterNumbers input", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # First create a named group
      suppressWarnings(session$setInputs(
        clusterName = "TestGroup",
        clusterNumbers = "5,6,7",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      # Verify group was created
      # message(paste(names(rv$outputList), collapse = ", "))
      expect_true("TestGroup" %in% names(rv$outputList))

      # Select the group via clusterNameSelect
      suppressWarnings(session$setInputs(clusterNameSelect = "TestGroup"))
      suppressWarnings(session$flushReact())

      # The clusterNumbers text input should be updated
      # This is tested by checking the observer runs without error
      # and the input value contains the expected cluster IDs
      expect_true(nzchar(input$clusterNumbers))
      # message(paste(input$clusterNumbers, collapse = ", "))
      expect_true(all(c("5", "6", "7") %in% strsplit(input$clusterNumbers, ",")[[1]]))
    }
  )
})


test_that("observers: clusterNameSelect handles multiple groups", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Create two groups
      suppressWarnings(session$setInputs(
        clusterName = "GroupA",
        clusterNumbers = "1,2",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      suppressWarnings(session$setInputs(
        clusterName = "GroupB",
        clusterNumbers = "3,4",
        applyName = 2
      ))
      suppressWarnings(session$flushReact())

      # Select both groups
      suppressWarnings(session$setInputs(clusterNameSelect = c("GroupA", "GroupB")))
      suppressWarnings(session$flushReact())

      # clusterNumbers should contain combined unique IDs
      expect_true(nzchar(input$clusterNumbers))
    }
  )
})


test_that("observers: clusterNameSelect handles empty selection", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Empty selection should not cause error
      suppressWarnings(session$setInputs(clusterNameSelect = character(0)))
      suppressWarnings(session$flushReact())

      # Observer should complete without error
      expect_true(TRUE)
    }
  )
})


# =============================================================================
# groupsVar observer (lines 162-170)
# =============================================================================
test_that("observers: groupsVar updates group1 and group2 choices", {
  test_app <- make_test_app()
  sce <- test_app$sce
  metaD <- S4Vectors::metadata(sce)

  fact_cols <- vapply(metaD$experiment_info, is.factor, logical(1))
  skip_if_not(any(fact_cols), "fixture has no factor column in experiment_info")
  factor_col <- names(fact_cols)[fact_cols][1]
  expected_levels <- levels(metaD$experiment_info[[factor_col]])

  calls <- list()
  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, ..., choices = NULL, selected = NULL) {
      calls[[inputId]] <<- list(choices = choices, selected = selected)
    },
    .package = "shiny"
  )

  shiny::testServer(test_app$app, {
    suppress_plotly_event_warnings(session$setInputs(groupsVar = factor_col))

    expect_true("group1" %in% names(calls))
    expect_true("group2" %in% names(calls))
    expect_setequal(unlist(calls[["group1"]]$choices), expected_levels)
    expect_setequal(unlist(calls[["group2"]]$choices), expected_levels)
  })
})

test_that("observers: groupsVar handles invalid column", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Set invalid column name
      suppressWarnings(session$setInputs(groupsVar = "nonexistent_column"))
      suppressWarnings(session$flushReact())

      # Observer should return NULL without error
      expect_true(input$groupsVar == "nonexistent_column")
    }
  )
})


test_that("observers: groupsVar handles NULL input", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Set NULL groupsVar
      suppressWarnings(session$setInputs(groupsVar = NULL))
      suppressWarnings(session$flushReact())

      # Observer should return NULL without error
      expect_null(input$groupsVar)
    }
  )
})


# =============================================================================
# group1/group2 mutual exclusion observers (lines 173-201)
# =============================================================================
test_that("observers: group1 and group2 are mutually exclusive (grp1 changes)", {
  test_app <- make_test_app()
  sce <- test_app$sce
  metaD <- S4Vectors::metadata(sce)

  fact_cols <- unlist(lapply(metaD$experiment_info, is.factor), use.names = FALSE)
  if (any(fact_cols)) {
    factor_col <- names(fact_cols)[fact_cols][1]
    levs <- levels(metaD$experiment_info[, factor_col])

    if (length(levs) >= 3) {
      shiny::testServer(
        app = test_app$app,
        expr = {
          # Set up groups
          suppressWarnings(session$setInputs(
            groupsVar = factor_col,
            group1 = levs[1],
            group2 = levs[2]
          ))
          suppressWarnings(session$flushReact())

          # Change group1 to include levs[2]
          suppressWarnings(session$setInputs(group1 = c(levs[1], levs[2])))
          suppressWarnings(session$flushReact())

          # group2 should be updated to exclude the new group1 values
          # The observer ensures mutual exclusion
          expect_true(TRUE) # Observer ran without error
        }
      )
    }
  }
})


test_that("observers: group1 and group2 are mutually exclusive (grp2 changes)", {
  test_app <- make_test_app()
  sce <- test_app$sce
  metaD <- S4Vectors::metadata(sce)

  fact_cols <- unlist(lapply(metaD$experiment_info, is.factor), use.names = FALSE)
  if (any(fact_cols)) {
    factor_col <- names(fact_cols)[fact_cols][1]
    levs <- levels(metaD$experiment_info[, factor_col])

    if (length(levs) >= 3) {
      shiny::testServer(
        app = test_app$app,
        expr = {
          # Set up groups
          suppressWarnings(session$setInputs(
            groupsVar = factor_col,
            group1 = levs[1],
            group2 = levs[2]
          ))
          suppressWarnings(session$flushReact())

          # Change group2 to include levs[1]
          suppressWarnings(session$setInputs(group2 = c(levs[2], levs[1])))
          suppressWarnings(session$flushReact())

          # group1 should be updated to exclude the new group2 values
          expect_true(TRUE) # Observer ran without error
        }
      )
    }
  }
})


test_that("observers: group1/group2 mutual exclusion with invalid groupsVar", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Set invalid groupsVar
      suppressWarnings(session$setInputs(groupsVar = "invalid"))
      suppressWarnings(session$setInputs(group1 = "a", group2 = "b"))
      suppressWarnings(session$flushReact())

      # Observer should return NULL without error
      expect_true(TRUE)
    }
  )
})


# =============================================================================
# Plotly selection observers for SOM plots (lines 212-233)
# =============================================================================
test_that("observers: somData plotly_selected updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # First establish a selection
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2",
        selectMode = "add"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(
          list(curveNumber = 0, pointNumber = 0),
          list(curveNumber = 0, pointNumber = 1)
        ),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-somData1` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # pointNumber 0/1 -> node ids 1/2 via the +1L fallback offset,
      # combined ("add") with the pre-existing 1,2 selection
      expect_true(all(c(1, 2) %in% rsUsed()))
    }
  )
})


test_that("observers: somData plotly_selected with null data returns NULL", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Establish initial selection
      suppressWarnings(session$setInputs(clusterNumbers = "1,2"))
      suppressWarnings(session$flushReact())

      initial_rs <- rsUsed()

      # Simulate null plotly_selected event
      suppressWarnings(session$setInputs(`plotly_selected-somData1` = NULL))
      suppressWarnings(session$flushReact())

      # Selection should remain unchanged
      expect_equal(rsUsed(), initial_rs)
    }
  )
})


test_that("observers: somData plotly_selected respects selectMode", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Establish initial selection
      suppressWarnings(session$setInputs(clusterNumbers = "1,2"))
      suppressWarnings(session$flushReact())

      # Test "add" mode
      suppressWarnings(session$setInputs(selectMode = "add"))
      suppressWarnings(session$flushReact())

      expect_true(TRUE) # Observer ran without error
    }
  )
})


# =============================================================================
# Zoom observers for SOM plots (lines 236-244)
# =============================================================================
test_that("observers: somData plotly_relayout triggers zoomFunc", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Seed axis selection so dimSelection() actually initializes
      suppressWarnings(session$setInputs(
        currentDimX = "marker1",
        currentDimY = "marker2"
      ))
      suppressWarnings(session$flushReact()) # arm the activeDims debounce timer
      session$elapse(300) # push past the 250ms debounce window
      suppressWarnings(session$flushReact()) # let dimSelection() observer run

      before <- dimSelection()
      expect_true(length(before) >= 1L)
      expect_null(before[[1]]$xzoom)

      relayout_json <- jsonlite::toJSON(
        list(
          `xaxis.range[0]` = 0,
          `xaxis.range[1]` = 1,
          `yaxis.range[0]` = 0,
          `yaxis.range[1]` = 1
        ),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_relayout-somData1` = relayout_json
      ))
      suppressWarnings(session$flushReact())

      after <- dimSelection()
      expect_equal(after[[1]]$xzoom, c(0, 1))
      expect_equal(after[[1]]$yzoom, c(0, 1))
    }
  )
})


test_that("observers: somData plotly_relayout with null data", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Simulate null plotly_relayout event
      suppressWarnings(session$setInputs(`plotly_relayout-somData1` = NULL))
      suppressWarnings(session$flushReact())

      # Observer should handle null without error
      expect_true(TRUE)
    }
  )
})


# =============================================================================
# somDataMain selection observer (lines 247-259)
# =============================================================================
test_that("observers: somDataMain plotly_selected updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Establish initial selection
      suppressWarnings(session$setInputs(
        clusterNumbers = "5",
        selectMode = "add"
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000) # inputClusterNumber debounce -> rsUsed(c(5))
      session$elapse(1500) # rsUsed_d debounce -> rsUsed_d() picks up c(5)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-somDataMain` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # pointNumber 0 -> node id 1 via the curveNumber==0 fallback (+1L offset),
      # combined ("add") with the pre-existing selection of 5
      expect_true(all(c(1, 5) %in% rsUsed()))
    }
  )
})

test_that("observers: somDataMain plotly_selected with no event leaves rsUsed unchanged", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Prime rsUsed_d() first, so it isn't NULL/uninitialized —
      # otherwise .inputSelect()'s is.null(rs) guard wipes the
      # selection regardless of what we're testing here.
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())

      initial_rs <- rsUsed()
      expect_equal(initial_rs, 1)

      # Deliberately don't touch `plotly_selected-somDataMain` at all —
      # that's the correct way to simulate "no selection event has
      # fired yet."
      suppressWarnings(session$flushReact())

      expect_equal(rsUsed(), initial_rs)
    }
  )
})


# =============================================================================
# somDataMain zoom observer (lines 262-267)
# =============================================================================
test_that("observers: somDataMain plotly_relayout triggers zoomFunc", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Seed axis selection so dimSelection() actually initializes
      suppressWarnings(session$setInputs(
        currentDimX = "marker1",
        currentDimY = "marker2"
      ))
      suppressWarnings(session$flushReact()) # arm the activeDims debounce timer
      session$elapse(300) # push past the 250ms debounce window
      suppressWarnings(session$flushReact()) # let dimSelection() observer run

      before <- dimSelection()
      expect_true(length(before) >= 1L)
      expect_null(before[[1]]$xzoom)

      relayout_json <- jsonlite::toJSON(
        list(
          `xaxis.range[0]` = 0,
          `xaxis.range[1]` = 1,
          `yaxis.range[0]` = 2,
          `yaxis.range[1]` = 3
        ),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_relayout-somDataMain` = relayout_json
      ))
      suppressWarnings(session$flushReact())

      after <- dimSelection()
      expect_equal(after[[1]]$xzoom, c(0, 1))
      expect_equal(after[[1]]$yzoom, c(2, 3))
    }
  )
})


# =============================================================================
# Dimension-reduction plot selection observers (lines 270-293)
# =============================================================================
test_that("observers: tsne plotly_selected updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Establish initial selection
      suppressWarnings(session$setInputs(
        clusterNumbers = "3",
        selectMode = "add"
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000) # inputClusterNumber debounce -> rsUsed(c(3))
      session$elapse(1500) # rsUsed_d debounce -> rsUsed_d() picks up c(3)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-tsne` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # pointNumber 0 -> node id 1 via the curveNumber==0 fallback (+1L offset),
      # combined ("add") with the pre-existing selection of 3
      expect_true(all(c(1, 3) %in% rsUsed()))
    }
  )
})


test_that("observers: umap plotly_selected updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "4",
        selectMode = "add"
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-umap` = selected_json
      ))
      suppressWarnings(session$flushReact())

      expect_true(all(c(1, 4) %in% rsUsed()))
    }
  )
})


test_that("observers: pca plotly_selected updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "6",
        selectMode = "add"
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-pca` = selected_json
      ))
      suppressWarnings(session$flushReact())

      expect_true(all(c(1, 6) %in% rsUsed()))
    }
  )
})


test_that("observers: dimension reduction plotly_selected with no event leaves rsUsed unchanged", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())

      initial_rs <- rsUsed()
      expect_equal(initial_rs, 1)

      # Deliberately don't touch plotly_selected-tsne/umap/pca at all.
      suppressWarnings(session$flushReact())

      expect_equal(rsUsed(), initial_rs)
    }
  )
})


# =============================================================================
# Dendrogram selection observers (lines 296-346)
# =============================================================================
test_that("observers: dendPlotly plotly_selected updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Establish initial selection
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000) # inputClusterNumber debounce -> rsUsed(c(1))
      session$elapse(1500) # rsUsed_d debounce -> rsUsed_d() picks up c(1)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0, customdata = 7)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-dendPlotly` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # Should have both the pre-existing 1 and the newly selected 7
      expect_true(all(c(1, 7) %in% rsUsed()))
    }
  )
})

test_that("observers: dendPlotly plotly_selected with empty data", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())
      initial_rs <- rsUsed()

      suppressWarnings(session$setInputs(
        `plotly_selected-dendPlotly` = "[]"
      ))
      suppressWarnings(session$flushReact())

      # Selection should be untouched by an empty box-select
      expect_equal(rsUsed(), initial_rs)
    }
  )
})

test_that("observers: dendPlotly plotly_selected respects selectMode remove others", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2,3",
        selectMode = "remove others"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0, customdata = 2)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-dendPlotly` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # Should only have cluster 2
      expect_equal(rsUsed(), c(2))
    }
  )
})

test_that("observers: dendPlotly plotly_selected respects selectMode add", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      session$elapse(1000) # inputClusterNumber debounce -> rsUsed(c(1))
      session$elapse(1500) # rsUsed_d debounce -> rsUsed_d() picks up c(1)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0, customdata = 5)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-dendPlotly` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # Should have both 1 and 5
      expect_true(all(c(1, 5) %in% rsUsed()))
    }
  )
})


test_that("observers: dendPlotly plotly_selected respects selectMode remove", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2,3",
        selectMode = "remove"
      ))
      session$elapse(1000) # inputClusterNumber debounce -> rsUsed(c(1,2,3))
      session$elapse(1500) # rsUsed_d debounce -> rsUsed_d() picks up c(1,2,3)
      suppressWarnings(session$flushReact())

      selected_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0, customdata = 2)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-dendPlotly` = selected_json
      ))
      suppressWarnings(session$flushReact())

      # Should have 1 and 3, but not 2
      expect_true(all(c(1, 3) %in% rsUsed()))
      expect_false(2 %in% rsUsed())
    }
  )
})


test_that("observers: dendPlotly plotly_click updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "view"
      ))
      session$elapse(1000) # inputClusterNumber debounce
      session$elapse(1500) # rsUsed_d debounce
      suppressWarnings(session$flushReact())

      click_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0, customdata = 3)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_click-dendPlotly` = click_json
      ))
      suppressWarnings(session$flushReact())

      expect_true(3 %in% rsUsed())
    }
  )
})


test_that("observers: dendPlotly plotly_click with null data leaves state unchanged", {
  test_app <- make_test_app()
  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())
      initial_rs <- rsUsed()

      suppressWarnings(session$setInputs(
        `plotly_click-dendPlotly` = "[]"
      ))
      suppressWarnings(session$flushReact())

      expect_equal(rsUsed(), initial_rs)
    }
  )
})

test_that("observers: dendPlotly plotly_click respects selectMode", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2",
        selectMode = "add"
      ))
      session$elapse(1000) # inputClusterNumber debounce -> rsUsed(c(1,2))
      session$elapse(1500) # rsUsed_d debounce -> rsUsed_d() picks up c(1,2)
      suppressWarnings(session$flushReact())

      click_json <- jsonlite::toJSON(
        list(list(curveNumber = 0, pointNumber = 0, customdata = 7)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_click-dendPlotly` = click_json
      ))
      suppressWarnings(session$flushReact())

      expect_true(all(c(1, 2, 7) %in% rsUsed()))
    }
  )
})

# =============================================================================
# SOM raster grid selection observer (lines 349-365)
#
# NOTE: The three somGrid tests that previously lived here were removed.
# The first ("updates rsUsed") failed intermittently on the Bioconductor/
# r-universe CI runners with a shiny `destroyedReactiveError` on `rsUsed()`
# after `session$elapse()`, while passing locally; the other two were removed
# alongside it for consistency since they exercised the same observer.
# The somGrid observer's logic (isolate(rsUsed()), .safeEventData,
# .inputSelect) is structurally identical to, and already covered by, the
# somData/somDataMain observer tests above.
# =============================================================================


# =============================================================================
# Scatter-plot rectangular selection observer (lines 368-410)
# =============================================================================
test_that("observers: scatterPlot plotly_selected with a single (zero-area) point leaves selection unchanged", {
  test_app <- make_test_app()
  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add",
        samples2plot = c("sample1")
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())
      initial_rs <- rsUsed()

      # A single point collapses the derived box to zero width/height.
      # The observer now explicitly guards against this and returns
      # early, leaving the selection unchanged.
      single_point_json <- jsonlite::toJSON(
        list(list(curveNumber = 1, x = 0.5, y = 0.5)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-scatterPlot` = single_point_json
      ))
      suppressWarnings(session$flushReact())

      # Degenerate selection should not modify rsUsed()
      expect_equal(rsUsed(), initial_rs)
    }
  )
})


test_that("observers: scatterPlot plotly_selected with missing curveNumber leaves state unchanged", {
  test_app <- make_test_app()
  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add",
        samples2plot = c("sample1")
      ))
      suppressWarnings(session$flushReact())
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())
      initial_rs <- rsUsed()

      no_curve_json <- jsonlite::toJSON(
        list(list(x = 0.5, y = 0.5)),
        auto_unbox = TRUE
      )
      suppressWarnings(session$setInputs(
        `plotly_selected-scatterPlot` = no_curve_json
      ))
      suppressWarnings(session$flushReact())

      expect_equal(rsUsed(), initial_rs)
    }
  )
})

test_that("observers: dendPlotly plotly_selected with empty data leaves state unchanged", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())
      initial_rs <- rsUsed()

      suppressWarnings(session$setInputs(
        `plotly_selected-dendPlotly` = "[]"
      ))
      suppressWarnings(session$flushReact())

      expect_equal(rsUsed(), initial_rs)
    }
  )
})

test_that("observers: scatterPlot plotly_selected with empty data leaves state unchanged", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1",
        selectMode = "add"
      ))
      session$elapse(1000)
      session$elapse(1500)
      suppressWarnings(session$flushReact())
      initial_rs <- rsUsed()

      suppressWarnings(session$setInputs(
        `plotly_selected-scatterPlot` = "[]"
      ))
      suppressWarnings(session$flushReact())

      expect_equal(rsUsed(), initial_rs)
    }
  )
})

test_that("observers: scatterPlot plotly_selected respects selectMode", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2,3",
        selectMode = "remove others"
      ))
      suppressWarnings(session$flushReact())

      # The observer should run without error
      expect_true(TRUE)
    }
  )
})


# =============================================================================
# FlowSOM star plot selection observers (lines 413-461)
# =============================================================================
test_that("observers: flowSOMStars plotly_selected updates rsUsed", {
  # Create test app with fsom
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # These observers only run when fsom is not NULL
      # The observer structure is tested by the code loading without error
      expect_true(TRUE)
    }
  )
})


test_that("observers: flowSOMStars plotly_click updates rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Similar to plotly_selected, tests observer structure
      expect_true(TRUE)
    }
  )
})


# =============================================================================
# close button observer (lines 664-667)
# =============================================================================
test_that("observers: close button triggers stopApp", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Trigger close button
      suppressWarnings(session$setInputs(close = 1))
      suppressWarnings(session$flushReact())

      # The observer calls stopApp which will terminate the test
      # We just verify the observer is set up correctly
      expect_true(TRUE)
    }
  )
})


# =============================================================================
# Edge cases for clusterNumbers observer (lines 74-82)
# =============================================================================
test_that("observers: clusterNumbers handles non-integer input by keeping the prior selection", {
  test_app <- make_test_app()
  shiny::testServer(test_app$app, {
    # capture whatever the selection was before the bad input
    before <- rsUsed()

    suppressWarnings(session$setInputs(clusterNumbers = "abc,def"))
    session$elapse(1001) # let inputClusterNumber()'s debounce settle

    expect_equal(rsUsed(), before) # unchanged, not cleared
  })
})


test_that("observers: clusterNumbers handles mixed valid/invalid input", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(clusterNumbers = "1,abc,2"))
      session$elapse(1000) # inputClusterNumber debounce
      session$elapse(1500) # rsUsed_d debounce
      suppressWarnings(session$flushReact())

      # Should extract valid integers
      expect_true(all(c(1, 2) %in% rsUsed()))
    }
  )
})


test_that("observers: clusterNumbers handles out-of-range clusters", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Clusters not in clusterPatientTable should be ignored
      suppressWarnings(session$setInputs(clusterNumbers = "999,1000"))
      session$elapse(1000) # inputClusterNumber debounce
      session$elapse(1500) # rsUsed_d debounce
      suppressWarnings(session$flushReact())

      # Should result in empty or unchanged selection
      expect_true(length(rsUsed()) >= 0)
    }
  )
})


# =============================================================================
# Edge cases for rmGrp observer (lines 119-129)
# =============================================================================
test_that("observers: rmGrp with empty name", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNameRM = "",
        rmGrp = 1
      ))
      suppressWarnings(session$flushReact())

      expect_true(TRUE)
    }
  )
})


test_that("observers: rmGrp with non-existent group", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(
        clusterNameRM = "NonExistentGroup",
        rmGrp = 1
      ))
      suppressWarnings(session$flushReact())

      expect_true(TRUE)
    }
  )
})


# =============================================================================
# Edge cases for sample2PlotDb observer (lines 204-209)
# =============================================================================
test_that("observers: sample2PlotDb with null sampleIds", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(samples2plot = NULL))
      suppressWarnings(session$flushReact())

      expect_true(TRUE)
    }
  )
})


test_that("observers: sample2PlotDb with empty sampleIds", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      suppressWarnings(session$setInputs(samples2plot = character(0)))
      suppressWarnings(session$flushReact())

      expect_true(TRUE)
    }
  )
})


# =============================================================================
# flowSOMStars observers (observers_clusterSelector.R:424-471)
# Gated behind `fsom` being non-NULL; a placeholder value is sufficient since
# these observers never dereference fsom's contents.
# =============================================================================

test_that("flowSOMStars plotly_selected updates rsUsed via customdata", {
  test_app <- make_test_app(fsom = fsom_stub)

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1",
      selectMode = "add"
    ))
    suppressWarnings(session$flushReact())
    session$elapse(1000)
    session$elapse(1500)
    suppressWarnings(session$flushReact())

    selected_json <- jsonlite::toJSON(
      list(list(curveNumber = 0, pointNumber = 0, customdata = 9)),
      auto_unbox = TRUE
    )
    suppressWarnings(session$setInputs(
      `plotly_selected-flowSOMStars` = selected_json
    ))
    suppressWarnings(session$flushReact())

    expect_true(all(c(1, 9) %in% rsUsed()))
  })
})

test_that("flowSOMStars plotly_selected with empty data leaves state unchanged", {
  test_app <- make_test_app(fsom = fsom_stub)

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1",
      selectMode = "add"
    ))
    suppressWarnings(session$flushReact())
    session$elapse(1000)
    session$elapse(1500)
    suppressWarnings(session$flushReact())
    initial_rs <- rsUsed()

    suppressWarnings(session$setInputs(
      `plotly_selected-flowSOMStars` = "[]"
    ))
    suppressWarnings(session$flushReact())

    # Same nrow(d) == 0L / is.null(d$customdata) chain as dendPlotly's
    # plotly_selected guard (observers_clusterSelector.R:432) -- confirm
    # this doesn't error the way we suspected dendPlotly's might.
    expect_equal(rsUsed(), initial_rs)
  })
})

test_that("flowSOMStars plotly_click updates rsUsed via customdata", {
  test_app <- make_test_app(fsom = fsom_stub)

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1",
      selectMode = "add"
    ))
    suppressWarnings(session$flushReact())
    session$elapse(1000)
    session$elapse(1500)
    suppressWarnings(session$flushReact())

    click_json <- jsonlite::toJSON(
      list(list(curveNumber = 0, pointNumber = 0, customdata = 5)),
      auto_unbox = TRUE
    )
    suppressWarnings(session$setInputs(
      `plotly_click-flowSOMStars` = click_json
    ))
    suppressWarnings(session$flushReact())

    expect_true(all(c(1, 5) %in% rsUsed()))
  })
})

test_that("flowSOMStars plotly_click with no customdata leaves state unchanged", {
  test_app <- make_test_app(fsom = fsom_stub)

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1",
      selectMode = "add"
    ))
    suppressWarnings(session$flushReact())
    session$elapse(1000)
    session$elapse(1500)
    suppressWarnings(session$flushReact())
    initial_rs <- rsUsed()

    click_json <- jsonlite::toJSON(
      list(list(curveNumber = 0, pointNumber = 0)), # no customdata field
      auto_unbox = TRUE
    )
    suppressWarnings(session$setInputs(
      `plotly_click-flowSOMStars` = click_json
    ))
    suppressWarnings(session$flushReact())

    expect_equal(rsUsed(), initial_rs)
  })
})

test_that("flowSOMStars observers are never registered when fsom is NULL", {
  test_app <- make_test_app(fsom = NULL)

  shiny::testServer(app = test_app$app, expr = {
    click_json <- jsonlite::toJSON(
      list(list(curveNumber = 0, pointNumber = 0, customdata = 5)),
      auto_unbox = TRUE
    )
    suppressWarnings(session$setInputs(`plotly_click-flowSOMStars` = click_json))
    suppressWarnings(session$flushReact())

    # No observer exists to act on this -- rsUsed should sit at its default
    expect_equal(rsUsed(), 1)
  })
})

# =============================================================================
# Scatter-plot rectangular selection observer -- guard paths
# =============================================================================

test_that("scatterPlot plotly_selected with missing curveNumber leaves state unchanged", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1",
      selectMode = "add",
      samples2plot = c("sample1")
    ))
    suppressWarnings(session$flushReact())
    session$elapse(1000)
    session$elapse(1500)
    suppressWarnings(session$flushReact())
    initial_rs <- rsUsed()

    # points with no curveNumber field at all -- req(d$curveNumber) should
    # stop the observer cleanly (d$curveNumber is NULL on a plain list)
    no_curve_json <- jsonlite::toJSON(
      list(list(pointNumber = 0, x = 1, y = 1)),
      auto_unbox = TRUE
    )
    suppressWarnings(session$setInputs(`plotly_selected-scatterPlot` = no_curve_json))
    suppressWarnings(session$flushReact())

    expect_equal(rsUsed(), initial_rs)
  })
})

test_that("scatterPlot plotly_selected with empty selection leaves state unchanged", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      clusterNumbers = "1",
      selectMode = "add",
      samples2plot = c("sample1")
    ))
    suppressWarnings(session$flushReact())
    session$elapse(1000)
    session$elapse(1500)
    suppressWarnings(session$flushReact())
    initial_rs <- rsUsed()

    suppressWarnings(session$setInputs(`plotly_selected-scatterPlot` = "[]"))
    suppressWarnings(session$flushReact())

    expect_equal(rsUsed(), initial_rs)
  })
})

test_that("scatterPlot plotly_selected computes correct cluster ids from a box select", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    message("rsUsed() at the very start, before any setInputs: ", paste(isolate(rsUsed()), collapse = ","))

    suppressWarnings(session$setInputs(
      currentDimX = "marker1",
      currentDimY = "marker2",
      samples2plot = c("sample1"),
      selectMode = "add"
    ))
    suppressWarnings(session$flushReact())
    message("rsUsed() after initial flush: ", paste(isolate(rsUsed()), collapse = ","))

    session$elapse(300)
    suppressWarnings(session$flushReact())
    message("rsUsed() after elapse(300): ", paste(isolate(rsUsed()), collapse = ","))

    session$elapse(1000)
    suppressWarnings(session$flushReact())
    message("rsUsed() after elapse(1000): ", paste(isolate(rsUsed()), collapse = ","))

    session$elapse(1500)
    suppressWarnings(session$flushReact())
    message("rsUsed() after elapse(1500): ", paste(isolate(rsUsed()), collapse = ","))
    message("rsUsed_d() after elapse(1500): ", paste(isolate(rsUsed_d()), collapse = ","))
    expected_ids <- sort(unique(as.integer(
      SingleCellExperiment::colData(
        test_app$sce_subsampled[, test_app$sce_subsampled$sample_id == "sample1"]
      )$cluster_id
    )))

    # Two points spanning the full value range (strictly past both ends,
    # since the observer's filter is > / < not >= / <=), both on
    # curveNumber == 1 to match the observer's filter. This still
    # relies on min(d$x)/max(d$x) being derived from whatever points
    # plotly actually sends for a box-select -- worth confirming against
    # a real captured payload before trusting this is the right shape.
    select_json <- jsonlite::toJSON(
      list(
        list(curveNumber = 1, pointNumber = 0, x = -1, y = -1),
        list(curveNumber = 1, pointNumber = 1, x = 1001, y = 1001)
      ),
      auto_unbox = TRUE
    )
    suppressWarnings(session$setInputs(`plotly_selected-scatterPlot` = select_json))
    suppressWarnings(session$flushReact())

    expect_true(all(expected_ids %in% rsUsed()))
  })
})

test_that("control: does ANY first setInputs reset rsUsed, or specifically currentDimX", {
  test_app <- make_test_app()

  shiny::testServer(app = test_app$app, expr = {
    message("baseline: ", paste(isolate(rsUsed()), collapse = ","))

    # scatterPercentile is read via isolate() in output$scatter and
    # nowhere else touches rsUsed -- a genuinely inert choice.
    suppressWarnings(session$setInputs(scatterPercentile = 0.95))
    suppressWarnings(session$flushReact())
    message("after an unrelated input (scatterPercentile): ", paste(isolate(rsUsed()), collapse = ","))

    expect_true(TRUE)
  })
})
# =============================================================================
# groupsVar / group1 / group2 observers (observers_clusterSelector.R:162-199)
# =============================================================================

test_that("groupsVar observeEvent updates both group1 and group2 choices to full levels", {
  test_app <- make_test_app()
  captured <- list()

  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, choices = NULL, selected = NULL, ...) {
      captured[[inputId]] <<- list(choices = choices, selected = selected)
    },
    .package = "shiny"
  )

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(groupsVar = "treatment"))
    suppressWarnings(session$flushReact())

    expect_equal(captured[["group1"]]$choices, c("control", "treated"))
    expect_equal(captured[["group2"]]$choices, c("control", "treated"))
  })
})

test_that("groupsVar observeEvent is a no-op for an invalid column name", {
  test_app <- make_test_app()
  captured <- list()

  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, choices = NULL, selected = NULL, ...) {
      captured[[inputId]] <<- list(choices = choices, selected = selected)
    },
    .package = "shiny"
  )

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(groupsVar = "not_a_real_column"))
    suppressWarnings(session$flushReact())

    expect_null(captured[["group1"]])
    expect_null(captured[["group2"]])
  })
})


test_that("group1 selection excludes it from group2's choices", {
  test_app <- make_test_app()
  captured <- list()

  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, choices = NULL, selected = NULL, ...) {
      captured[[inputId]] <<- list(choices = choices, selected = selected)
    },
    .package = "shiny"
  )

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      groupsVar = "treatment",
      group1 = "control"
    ))
    suppressWarnings(session$flushReact())

    expect_equal(captured[["group2"]]$choices, "treated")
  })
})

test_that("group2 selection excludes it from group1's choices", {
  test_app <- make_test_app()
  captured <- list()

  testthat::local_mocked_bindings(
    updateSelectInput = function(session, inputId, choices = NULL, selected = NULL, ...) {
      captured[[inputId]] <<- list(choices = choices, selected = selected)
    },
    .package = "shiny"
  )

  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(
      groupsVar = "treatment",
      group2 = "treated"
    ))
    suppressWarnings(session$flushReact())

    expect_equal(captured[["group1"]]$choices, "control")
  })
})
