# CySA Shiny App Optimization Plan

Based on `profile1.Rprofvis` (993 980 samples, 20 ms interval, ~5.5 h total session).

## Context

- **99.1 % of time is unattributed** to package source (`promises::with_promise_domain`, `domain$wrapSync`, `withCallingHandlers`, `force`, etc.). This is mostly Shiny async/event-loop overhead and idle time.
- **0.9 % (~183 s) is file-specific work** and is the actionable portion.
- Dominant file-specific costs are plotting, rendering, and reactive recomputation.

Goals are therefore:
1. Reduce the amount of work done per render.
2. Reduce reactive invalidations and redundant recomputation.
3. Improve Shiny dependency tracking so less time is spent in the event loop.
4. Remove the archived `ggthemr` dependency.

## Success metrics

- Reduce combined `somPlot` + `ggplotly` + `raster::plot` + `points` time by at least 50 % in a fresh `profvis` run over the same workflow.
- Reduce `dimSelection` observer cost to negligible (< 0.1 s total).
- Reduce `tsne`/`umap` recomputation to zero for identical inputs.
- `devtools::check()` passes with no new warnings or errors.
- Smoke test (`tests/testthat/test-smoke.R`) still passes.

## Validation tools

- `profvis::profvis()` around the same user workflow.
- `options(shiny.reactlog = TRUE)` + `reactlog::reactlog_show()` for reactive dependency tracing.
- `system.time()` micro-benchmarks on individual functions.
- `devtools::load_all()` + `devtools::test()` + `devtools::check()`.

---

## Phase 1 — Quick wins (low risk, high impact)

### 1.1 Pre-compute axis limits instead of scanning `assays(sce)[[1]]`

**Problem:** `R/clusterSelector.R:1909-1912` computes `min`/`max` of `assays(sce)[[1]]` for every plot whenever axes change. On a large assay this is expensive and happens repeatedly.

**Change:**
- In `prepClusterSelectorData()` compute a limits data frame/matrix `channel_limits` with `min` and `max` per channel used in `dList`.
- Pass `channel_limits` into `clusterSelector()`.
- Replace the `dimSelection` observer body with a lookup:

```r
# in prepClusterSelectorData or clusterSelector start
channel_limits <- lapply(colsUsed, function(ch) {
  c(min = min(SOM_codes[, ch]), max = max(SOM_codes[, ch]))
})
names(channel_limits) <- colsUsed

# in dimSelection observer
dimSelection[[idx]] <- list(
  dims = c(d1, d2),
  xlim = c(channel_limits[[d1]]["min"], channel_limits[[d1]]["max"]),
  ylim = c(channel_limits[[d2]]["min"], channel_limits[[d2]]["max"]),
  xzoom = c(NULL, NULL),
  yzoom = c(NULL, NULL)
)
```

**Affected files:** `R/prepClusterSelector.R`, `R/clusterSelector.R`.

**Validation:**
- `profvis` shows no `assays` samples in the `dimSelection` observer.
- Micro-benchmark: `bench::mark()` the old vs new observer; new should be < 5 ms per plot.

### 1.2 Remove `ggthemr` / `set_swatch` from `drawProjection()`

**Problem:** `R/shinyAppFunctions.R:300` calls `ggthemr('flat')` and `R/shinyAppFunctions.R:316` calls `set_swatch()` on every render. `ggthemr` is archived on CRAN and slow.

**Change:**
- Replace `ggthemr('flat')` with explicit `theme_minimal()` or `theme_bw()`.
- Replace `set_swatch()` with a manual `scale_color_manual(values = mycolors)`.
- Remove `ggthemr` from `DESCRIPTION` `Imports` after all calls are gone.

**Affected files:** `R/shinyAppFunctions.R`, `DESCRIPTION`, `NAMESPACE` (regenerate with `devtools::document()`).

**Validation:**
- `profvis` shows no `ggthemr` or `set_swatch` samples.
- Visual regression: launch the app and confirm group colors look the same.
- `devtools::check()` no longer warns about `ggthemr`.

### 1.3 Share a single base `plotSOMScatter` object across `output$somData1-6`

**Problem:** `R/clusterSelector.R:1343-1433` repeats the same `plotSOMScatter()` + `somPlot()` + `ggplotly()` pattern six times. Each render builds the base plot from scratch even though only the channel pair and limits differ.

**Change:**
- Create a reactive that builds the common base data/ggplot once per `rs`:

```r
somBaseData <- reactive({
  rs <- rsUsed_d()
  req(rs)
  list(
    data = S4Vectors::metadata(sce)[[somCodesName]],
    stats = S4Vectors::metadata(sce)$SOM_stats,
    rs = rs
  )
}) %>% bindCache(rsUsed_d())
```

- Add a helper `makeSomPlot(baseData, dimPair, colorbyGroups, showGroups)`.
- Replace the six `renderPlotly` blocks with a loop or a single parameterized output.

**Affected files:** `R/clusterSelector.R`, `R/shinyAppFunctions.R`.

**Validation:**
- `profvis` shows only one `plotSOMScatter` sample per `rs` change, not six.
- All six plots still render and respond to selections.

### 1.4 Cache dimensionality reductions with `bindCache()`

**Problem:** `R/clusterSelector.R:1452-1474` defines `tsne`, `umap`, and `pca` reactives with `debounce()` but no cache. Unrelated invalidations cause recomputation.

**Change:**

```r
tsne <- reactive({
  tsneFunc(dimRedSelection = dimRedSelection(), perplexity = input$perplexity, sce, somCodesName)
}) %>%
  debounce(1000) %>%
  bindCache(dimRedSelection(), input$perplexity)
```

Do the same for `umap` and `pca` with their respective inputs.

**Affected files:** `R/clusterSelector.R`.

**Validation:**
- Toggle an unrelated input (e.g. `showlegend`) and verify `tsne`/`umap`/`pca` do not recompute (via `shiny.reactlog`).
- `profvis` shows only one `Rtsne::Rtsne` and one `umap::umap` sample per unique parameter set.

### 1.5 Return ggplot objects directly from `renderPlot()` where possible

**Problem:** `R/clusterSelector.R:1724`, `1784`, `1816` wrap plots in `print()`. This adds overhead and can force extra grid evaluation.

**Change:**
- Replace `renderPlot({ print(vlnPlot()) })` with `renderPlot({ vlnPlot() })`, etc.

**Affected files:** `R/clusterSelector.R`.

**Validation:**
- Plots still render.
- Fewer `print` samples in `profvis`.

---

## Phase 2 — Reactivity refactor (medium risk, higher impact)

### 2.1 Convert `outputList` from external environment to `reactiveValues()`

**Problem:** `outputList` is stored in an external `env` and accessed via `get()`/`assign()` in many hot paths (`R/clusterSelector.R:763,787,797,735,etc.`). This breaks Shiny's reactive dependency tracking and forces manual invalidation.

**Change:**
- In the server function create `rv <- reactiveValues(outputList = outputList)`.
- Replace all `get("outputList", envir = env)` with `rv$outputList`.
- Replace all `assign(x = "outputList", value = ..., envir = env)` with `rv$outputList <- ...`.
- Update `updatedoutputList()` to read `rv$outputList` and trigger off it.
- Keep the `env` argument for backward compatibility, but synchronize the value back to `env` on app close if callers still expect it.

**Affected files:** `R/clusterSelector.R`.

**Validation:**
- App launches and group creation/removal still works.
- `reactlog` shows proper reactive edges for `outputList`.
- `profvis` shows fewer `force` / `withCallingHandlers` cycles tied to manual invalidation.

### 2.2 Cache `drawProjection()` per view

**Problem:** `R/shinyAppFunctions.R:293` recomputes the group-colored projection for every SOM plot when `showGroups` is TRUE, including `ggplot_build`, `left_join`, and color scale setup.

**Change:**
- Pre-compute a data frame containing all SOM nodes, their x/y coordinates per channel pair, and group membership for the current `outputList`.
- Cache it with `bindCache(rsUsed_d(), colorbyGroups, showGroups, dimSelection())`.
- `drawProjection()` becomes a simple ggplot construction from the cached data frame.

**Affected files:** `R/shinyAppFunctions.R`, `R/clusterSelector.R`.

**Validation:**
- Switch between `somData1` and `somData6` without changing `rs`; `drawProjection` should not recompute.
- `profvis` shows one `drawProjection` sample per `rs`/`colorbyGroups` change, not per plot.

### 2.3 Defer heavy plots until their boxes are opened

**Problem:** `VlnPlot`, `VlnPlot2`, and `UpSet` reactives compute on startup even though their boxes start collapsed.

**Change:**
- Add an input or JavaScript observer that tracks which `shinydashboardPlus::box` is open.
- Wrap heavy reactives with `req(boxOpen("violin"))` so they only compute when visible.
- Alternative: use `shinyjs::hidden()`/`toggle()` or observe the `collapsed` state via `input$<boxId>_status` if the box has an id.

**Affected files:** `R/clusterSelector.R` (UI and server).

**Validation:**
- App startup time is faster.
- `profvis` taken immediately after launch shows no `vlnPlot`/`upSetPlot` samples until the box is expanded.

---

## Phase 3 — Rendering overhaul (higher risk, bigger payoff)

### 3.1 Replace `raster::plot()` + `points()` with a lighter visualization

**Problem:** `R/clusterSelector.R:1680` calls `raster::plot(res[[1]], addfun = res[[2]], maxnl=80)` and `R/clusterSelector.R:1713-1714` draws two `points()` overlays. This is the single most expensive user operation (22.8 s + 19.9 s).

**Change options:**

**Option A (recommended first step):** Pre-render the static raster background once as a base64 PNG or `annotation_raster()` layer, then overlay selected nodes with `geom_point()`.

**Option B:** Convert the raster to a data frame and plot with `ggplot2::geom_raster()` + `geom_point()` for selections, avoiding base graphics entirely.

**Option C:** Keep base raster but render it only once per session and invalidate only the overlay.

**Affected files:** `R/clusterSelector.R`.

**Validation:**
- `profvis` shows `raster::plot` + `points` time reduced by at least 60 %.
- The SOM raster still displays and selected nodes are highlighted.
- Compare screenshots before and after.

### 3.2 Vectorize violin data construction

**Problem:** `R/shinyAppFunctions.R:241-289` and `R/clusterSelector.R:1726-1781` build violin data by looping over `upsetSelection` and `rbind()`ing rows.

**Change:**
- Pre-subset `sce@metadata[[somCodesName]]` to selected markers once.
- Build a long data frame with `data.table::rbindlist()` or `tidyr::pivot_longer()` instead of the loop + `rbind()`.
- Consider computing violin data once per `upsetSelection` + `violinSelection` and caching it.

**Affected files:** `R/shinyAppFunctions.R`, `R/clusterSelector.R`.

**Validation:**
- `profvis` shows no `rbind` samples inside violin construction.
- Violin plots render identically for the same selections.

---

## Validation protocol for each phase

1. **Baseline:** Run the app through a representative workflow and capture a new `profvis` file. Record total file-specific time and the top 10 hotspots.
2. **Implement one phase at a time.**
3. **Test:**
   - `devtools::load_all()`
   - `devtools::test()`
   - `devtools::check()`
4. **Profile again** under the same workflow.
5. **Compare:** Report the change in file-specific time and in the targeted hotspots.

## Recommended order of work

1. Phase 1.1 (axis limits) — removes a recurring 1 s+ cost with minimal code change.
2. Phase 1.2 (remove `ggthemr`) — reduces render time and future-proofs dependencies.
3. Phase 1.5 (remove `print()` wrappers) — trivial.
4. Phase 1.4 (`bindCache` for DR) — removes repeated t-SNE/UMAP cost.
5. Phase 1.3 (shared base SOM plot) — biggest reduction in repeated plotting.
6. Phase 2.1 (`reactiveValues` for `outputList`) — improves reactivity hygiene.
7. Phase 2.2 (cache `drawProjection`) — complements 1.3.
8. Phase 2.3 (defer heavy plots) — improves startup and reduces idle invalidations.
9. Phase 3.1 (raster rendering) — addresses the single largest hotspot.
10. Phase 3.2 (vectorize violins) — lower priority until raster work is done.

## Notes

- `NAMESPACE` is generated by roxygen2. After changing imports, run `devtools::document()` and do not edit `NAMESPACE` by hand.
- `README.md` is generated from `README.Rmd`; edit `README.Rmd` and run `devtools::build_readme()` if the plan changes user-facing behavior.
- Keep the `env` argument and `outputList` mutation contract for backward compatibility until callers are updated.
