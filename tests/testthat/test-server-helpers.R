# Tests for server-helpers.R
# These test the non-reactive helper functions used in the clusterSelector server.


# =============================================================================
# .updateOutputListInputs()
# =============================================================================

test_that(".updateOutputListInputs() runs without error", {
  # Create mock session using shiny's internal mock session
  input <- shiny::reactiveValues(
    clusterNameSelect = NULL,
    compareStatsTo = NULL,
    relativeTo = NULL,
    upsetSelection = NULL,
    colorbyGroups = NULL,
    groupRM = NULL
  )

  outputList <- list(Rest = c(1L, 2L, 3L), selected = c(1L, 2L))

  # Create minimal metadata
  metaD <- list(
    experiment_info = data.frame(
      sample_id = c("S1", "S2"),
      age = c(25, 30),
      stringsAsFactors = FALSE
    )
  )

  # Use testServer to get a valid session
  shiny::testServer(
    app = shiny::shinyApp(ui = shiny::fluidPage(), server = function(input, output, session) {
      # Call the function with the real session
      expect_no_error({
        .updateOutputListInputs(session, input, outputList, metaD)
      })
    }),
    expr = {}
  )
})


# =============================================================================
# .computeRelativeCounts()
# =============================================================================

test_that(".computeRelativeCounts() with relativeToCol='none' returns percentages", {
  clusterPatientTable <- matrix(
    c(100, 200, 150, 50, 100, 100),
    nrow = 3, ncol = 2,
    dimnames = list(c("Sample1", "Sample2", "Sample3"), c("1", "2"))
  )

  rs <- c("1")

  expInfo <- data.frame(
    sample_id = c("Sample1", "Sample2", "Sample3"),
    stringsAsFactors = FALSE
  )

  outputList <- list()

  result <- .computeRelativeCounts(clusterPatientTable, rs, "none", expInfo, outputList)

  # Should be percentage of selected cells relative to total
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(result >= 0 & result <= 100))
})


test_that(".computeRelativeCounts() with numeric column uses expInfo", {
  clusterPatientTable <- matrix(
    c(100, 200, 150, 50, 100, 100),
    nrow = 3, ncol = 2,
    dimnames = list(c("S1", "S2", "S3"), c("1", "2"))
  )

  rs <- c("1")

  expInfo <- data.frame(
    total_cells = c(200, 400, 300),
    row.names = c("S1", "S2", "S3")
  )

  outputList <- list()

  result <- .computeRelativeCounts(clusterPatientTable, rs, "total_cells", expInfo, outputList)

  expect_type(result, "double")
  expect_length(result, 3)
  # 100/200*100 = 50, 200/400*100 = 50, 150/300*100 = 50
  expect_equal(unname(result), c(50, 50, 50))
})


test_that(".computeRelativeCounts() with outputList group uses group counts", {
  clusterPatientTable <- matrix(
    c(100, 200, 150, 50, 100, 100),
    nrow = 3, ncol = 2,
    dimnames = list(c("S1", "S2", "S3"), c("1", "2"))
  )

  rs <- c("1")

  expInfo <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    stringsAsFactors = FALSE
  )

  outputList <- list(RefGroup = c("1", "2")) # Reference includes both clusters

  result <- .computeRelativeCounts(clusterPatientTable, rs, "RefGroup", expInfo, outputList)

  expect_type(result, "double")
  expect_length(result, 3)
  # 100/150*100, 200/300*100, 150/250*100
  expect_equal(unname(result), c(100 / 150 * 100, 200 / 300 * 100, 150 / 250 * 100))
})


test_that(".computeRelativeCounts() ignores stale node ids in a saved group", {
  # Regression: a saved group referencing a node id that no longer exists
  # (e.g. saved before revalidation, or from a different SOM) crashed with
  # "subscript out of bounds" instead of being ignored.
  clusterPatientTable <- matrix(
    c(100, 200, 150, 50, 100, 100),
    nrow = 3, ncol = 2,
    dimnames = list(c("S1", "S2", "S3"), c("1", "2"))
  )
  rs <- c("1")
  expInfo <- data.frame(sample_id = c("S1", "S2", "S3"), stringsAsFactors = FALSE)
  outputList <- list(RefGroup = c("1", "2", "99999"))

  result <- expect_no_error(
    .computeRelativeCounts(clusterPatientTable, rs, "RefGroup", expInfo, outputList)
  )
  # Same result as if "99999" had never been in the group.
  expect_equal(unname(result), c(100 / 150 * 100, 200 / 300 * 100, 150 / 250 * 100))
})


test_that(".computeRelativeCounts() ignores a stale id in rs itself", {
  clusterPatientTable <- matrix(
    c(100, 200, 150, 50, 100, 100),
    nrow = 3, ncol = 2,
    dimnames = list(c("S1", "S2", "S3"), c("1", "2"))
  )
  expInfo <- data.frame(sample_id = c("S1", "S2", "S3"), stringsAsFactors = FALSE)

  result <- expect_no_error(
    .computeRelativeCounts(clusterPatientTable, c("1", "99999"), "none", expInfo, list())
  )
  expect_equal(unname(result), c(100 / 150 * 100, 200 / 300 * 100, 150 / 250 * 100))
})


# =============================================================================
# .validNodeCols()
# =============================================================================

test_that(".validNodeCols() keeps only ids that are columns of clusterPatientTable", {
  cpt <- matrix(1, nrow = 1, ncol = 3, dimnames = list("S1", c("1", "2", "5")))
  expect_equal(.validNodeCols(c(1L, 2L, 5L, 99L), cpt), c(1L, 2L, 5L))
  expect_equal(.validNodeCols(c(99L, NA_integer_), cpt), integer(0))
  expect_equal(.validNodeCols(integer(0), cpt), integer(0))
})


test_that(".computeRelativeCounts() errors when relativeToCol not found", {
  clusterPatientTable <- matrix(
    c(100, 200, 150, 50, 100, 100),
    nrow = 3, ncol = 2,
    dimnames = list(c("S1", "S2", "S3"), c("1", "2"))
  )

  rs <- c("1")
  expInfo <- data.frame(sample_id = c("S1", "S2", "S3"))
  outputList <- list()

  expect_error(
    {
      .computeRelativeCounts(clusterPatientTable, rs, "NonExistent", expInfo, outputList)
    },
    "relativeToCol 'NonExistent' not found"
  )
})


# =============================================================================
# .rebuildOutputList()
# =============================================================================

test_that(".rebuildOutputList() rebuilds Rest group correctly", {
  outputList <- list(
    Group1 = c(1L, 2L),
    Group2 = c(3L, 4L)
  )

  sceLevels <- as.character(1:10)

  result <- .rebuildOutputList(outputList, sceLevels)

  expect_equal(result$Group1, c(1L, 2L))
  expect_equal(result$Group2, c(3L, 4L))
  expect_equal(result$Rest, 5:10)
})


test_that(".rebuildOutputList() removes empty groups", {
  outputList <- list(
    Group1 = c(1L, 2L),
    EmptyGroup = integer(0)
  )

  sceLevels <- as.character(1:5)

  result <- .rebuildOutputList(outputList, sceLevels)

  expect_true("Group1" %in% names(result))
  expect_false("EmptyGroup" %in% names(result))
})


test_that(".rebuildOutputList() with removed clusters updates all groups", {
  outputList <- list(
    Group1 = c(1L, 2L, 5L),
    Group2 = c(3L, 4L, 5L)
  )

  sceLevels <- as.character(1:6)
  removed <- c(5L, 6L)

  result <- .rebuildOutputList(outputList, sceLevels, removed)

  expect_equal(result$Group1, c(1L, 2L))
  expect_equal(result$Group2, c(3L, 4L))
  # Rest should be clusters not in any group after removal
  # After removing 5,6 from groups: Group1=c(1,2), Group2=c(3,4)
  # Rest = sceLevels not in (1,2,3,4) = c(5,6)
  expect_equal(result$Rest, c(5L, 6L))
})


test_that(".rebuildOutputList() handles integer sceLevels", {
  outputList <- list(Group1 = c(1L, 2L))
  sceLevels <- 1:5

  result <- .rebuildOutputList(outputList, sceLevels)

  expect_equal(result$Group1, c(1L, 2L))
  expect_equal(result$Rest, 3:5)
})


# =============================================================================
# .inputSelect()
# =============================================================================

test_that(".inputSelect() with mode='view' returns node_ids", {
  d <- data.frame(key = c(1, 2, 3))
  rs <- c(4L, 5L)

  result <- .inputSelect(d, rs, "view")

  expect_equal(result, c(1L, 2L, 3L))
})


test_that(".inputSelect() with mode='add' combines rs and node_ids", {
  d <- data.frame(key = c(1, 2))
  rs <- c(3L, 4L)

  result <- .inputSelect(d, rs, "add")

  expect_equal(sort(result), c(1L, 2L, 3L, 4L))
})


test_that(".inputSelect() with mode='remove' subtracts node_ids from rs", {
  d <- data.frame(key = c(1, 2))
  rs <- c(1L, 2L, 3L, 4L)

  result <- .inputSelect(d, rs, "remove")

  expect_equal(sort(result), c(3L, 4L))
})


test_that(".inputSelect() with mode='remove others' returns intersection", {
  d <- data.frame(key = c(1, 2, 5))
  rs <- c(2L, 3L, 4L)

  result <- .inputSelect(d, rs, "remove others")

  expect_equal(result, c(2L))
})


test_that(".inputSelect() with NULL d returns rs", {
  result <- .inputSelect(NULL, c(1L, 2L), "view")
  expect_equal(result, c(1L, 2L))
})


test_that(".inputSelect() with NULL rs returns integer(0)", {
  d <- data.frame(key = c(1, 2))
  result <- .inputSelect(d, NULL, "view")
  expect_equal(result, integer(0))
})


test_that(".inputSelect() uses customdata when key not available", {
  d <- data.frame(customdata = c(10, 20, 30))
  rs <- c(1L)

  result <- .inputSelect(d, rs, "view")

  expect_equal(result, c(10L, 20L, 30L))
})


test_that(".inputSelect() uses pointNumber as last resort", {
  d <- data.frame(pointNumber = c(0, 1, 2), curveNumber = c(0, 0, 0))
  rs <- c(1L)

  result <- .inputSelect(d, rs, "view")

  # pointNumber is 0-indexed, so +1
  expect_equal(result, c(1L, 2L, 3L))
})


test_that(".inputSelect() filters to curveNumber == 0 when using pointNumber", {
  d <- data.frame(
    pointNumber = c(0, 1, 5, 6),
    curveNumber = c(0, 0, 1, 1)
  )
  rs <- c(1L)

  result <- .inputSelect(d, rs, "view")

  # Only curveNumber == 0 should be used
  expect_equal(result, c(1L, 2L))
})


test_that(".inputSelect() handles invalid node_ids gracefully", {
  d <- data.frame(key = c(NA, -1, 0, 5))
  rs <- c(1L)

  result <- .inputSelect(d, rs, "view")

  # NA, -1, and 0 should be filtered out
  expect_equal(result, c(5L))
})


test_that(".inputSelect() with all invalid node_ids returns rs", {
  d <- data.frame(key = c(NA, -1, 0))
  rs <- c(1L, 2L)

  result <- .inputSelect(d, rs, "view")

  expect_equal(result, c(1L, 2L))
})

test_that(".rebuildOutputList() recomputes a stale Rest instead of dropping it", {
    # Regression: the previous Rest counted as "used", so after naming a
    # group Rest became empty and disappeared from all group pickers.
    ol <- list(Rest = 1:10, selected = 3:4, grpA = 3:4)
    result <- .rebuildOutputList(ol, as.character(1:10))
    expect_equal(result$Rest, c(1L, 2L, 5:10))
    expect_equal(result$grpA, 3:4)
})
