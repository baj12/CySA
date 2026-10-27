test_that("server initialises reactive outputList with Rest group", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      expect_true("Rest" %in% names(rv$outputList))
      expect_length(rv$outputList$Rest, 10)
    }
  )
})

test_that("server updates rv$outputList when group is named", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Initial state: only Rest.
      expect_equal(sort(names(rv$outputList)), "Rest")

      # Simulate typing a name and applying a cluster selection.
      suppressWarnings(session$setInputs(
        clusterName = "MyGroup",
        clusterNumbers = "1,2",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      expect_true("MyGroup" %in% names(rv$outputList))
      expect_true(all(rv$outputList$MyGroup == c(1L, 2L)))
    }
  )
})

test_that("server reactive values are isolated from external env changes", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      expect_s4_class(test_app$sce, "SingleCellExperiment")
      expect_true("Rest" %in% names(rv$outputList))
      expect_type(dListRV(), "list")
    }
  )
})

test_that("server removes a single group when rmGrp is triggered", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # First, create a named group
      suppressWarnings(session$setInputs(
        clusterName = "TestGroup",
        clusterNumbers = "1,2",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      # Verify group was created
      expect_true("TestGroup" %in% names(rv$outputList))

      # Set the group to remove (clusterNameRM) and trigger rmGrp
      suppressWarnings(session$setInputs(
        clusterNameRM = "TestGroup",
        rmGrp = 1 # Increment to trigger the observer
      ))
      suppressWarnings(session$flushReact())

      # Verify group was removed
      expect_false("TestGroup" %in% names(rv$outputList))
    }
  )
})

test_that("server rmGrp does nothing when no group is selected", {
  test_app <- make_test_app()
  suppress_plotly_event_warnings({
    shiny::testServer(
      app = test_app$app,
      expr = {
        # Initial state: only Rest
        session$flushReact()

        initial_names <- names(rv$outputList)

        # Trigger rmGrp without selecting a group to remove
        suppressWarnings(session$setInputs(rmGrp = 1))
        suppressWarnings(session$flushReact())

        # State should be unchanged
        expect_equal(names(rv$outputList), initial_names)
      }
    )
  })
})

test_that("server rmGrp does nothing when selected group doesn't exist", {
  test_app <- make_test_app()

  suppress_plotly_event_warnings({
    shiny::testServer(
      app = test_app$app,
      expr = {
        session$flushReact()
        # Initial state: only Rest
        initial_names <- names(rv$outputList)

        # Try to remove a non-existent group
        suppressWarnings(session$setInputs(
          clusterNameRM = "NonExistentGroup",
          rmGrp = 1
        ))
        suppressWarnings(session$flushReact())

        # State should be unchanged
        expect_equal(names(rv$outputList), initial_names)
      }
    )
  })
})

# node-properties: reactives, renderer and selection ----

test_that("nodePropsDf combines stats, labels and ambiguity", {
  fx <- make_labeled_test_app()
  suppress_plotly_event_warnings({
    quiet_plotly_test(shiny::testServer(
      app = fx$app,
      expr = {
        session$flushReact()
        df <- st$nodePropsDf()
        expect_true(all(c("id", "n", "purity", "cell_type_majority",
          "medianRatio", "node_dispersion") %in% colnames(df)))
        expect_identical(nrow(df), 10L)
      }
    ))
  })
})

test_that("nodePropsPlot renders for labeled and unlabeled apps", {
  fx <- make_labeled_test_app()
  suppress_plotly_event_warnings({
    quiet_plotly_test(shiny::testServer(
      app = fx$app,
      expr = {
        session$flushReact()
        session$setInputs(nodePropsX = "purity", nodePropsY = "n")
        out <- output$nodePropsPlot
        expect_true(!is.null(out))
      }
    ))
  })
  fx2 <- make_test_app()
  suppress_plotly_event_warnings({
    quiet_plotly_test(shiny::testServer(
      app = fx2$app,
      expr = {
        session$flushReact()
        session$setInputs(nodePropsX = "medianRatio", nodePropsY = "n")
        out <- output$nodePropsPlot
        expect_true(!is.null(out))
      }
    ))
  })
})
