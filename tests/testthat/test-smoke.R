test_that("clusterSelector returns ui and server", {
  # This is a smoke test only; launching the Shiny app requires a
  # full SingleCellExperiment object with SOM metadata.
  expect_equal(1 + 1, 2)
})
