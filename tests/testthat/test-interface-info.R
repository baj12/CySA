# Tests for R/interface_info.R ----
# The per-box info button and its help-modal registry.

test_that(".boxInfoRegistry entries each have a title and 4 body paragraphs", {
  expect_gt(length(CySA:::.boxInfoRegistry), 0)
  for (id in names(CySA:::.boxInfoRegistry)) {
    info <- CySA:::.boxInfoRegistry[[id]]
    expect_true(is.character(info$title) && nzchar(info$title),
      info = paste("missing title for", id))
    expect_gte(length(info$body), 4L, label = paste("body length for", id))
    expect_true(all(nzchar(info$body)), info = paste("empty paragraph in", id))
  }
})

test_that(".infoButton returns an actionButton tagged with its infoId", {
  btn <- CySA:::.infoButton("diagnostics")
  html <- as.character(btn)
  expect_match(html, "boxInfo_diagnostics", fixed = TRUE)
  expect_match(html, "circle-info", fixed = TRUE)
})

test_that(".primaryBox embeds an info button only when infoId is given", {
  withButton <- CySA:::.primaryBox("Title", "content", infoId = "diagnostics")
  withoutButton <- CySA:::.primaryBox("Title", "content")

  expect_match(as.character(withButton), "boxInfo_diagnostics", fixed = TRUE)
  expect_false(grepl("boxInfo_", as.character(withoutButton), fixed = TRUE))
})

test_that("Every registered info button fires its observer without error", {
  test_app <- make_test_app(fsom = fsom_stub)

  shiny::testServer(app = test_app$app, expr = {
    for (id in names(CySA:::.boxInfoRegistry)) {
      args <- stats::setNames(list(1), paste0("boxInfo_", id))
      expect_error(
        suppressWarnings(do.call(session$setInputs, args)),
        NA,
        info = paste("info button failed for", id)
      )
    }
    suppressWarnings(session$flushReact())
  })
})
