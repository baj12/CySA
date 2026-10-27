# CySA 0.99.41

## New features

- Known cell types can now feed the app end to end. `prepClusterSelectorData()`
  and `prepClusterSelector()` accept a new `cellTypeCol` argument naming a
  colData column of per-cell labels; when given, each SOM node is scored
  against those labels: `metadata(sce)$node_labels` gains a data frame keyed
  by node with the node's majority label (`cell_type_majority`) and the share
  of its labeled cells carrying that label (`purity`), with `NA` for nodes
  holding no labeled cells. The label column is also carried onto the
  subsample. Independently of labels, `metadata(sce)$SOM_stats` gains a
  `node_dispersion` column -- the mean within-node marker SD, in the same
  cost class as the existing per-node stats.
- New `nodeAmbiguity()` computes, per node on the subsample, how close its
  cells sit to rival nodes (`medianRatio`, `pctAmbiguous`) -- the
  label-free counterpart of a mixed node.
- New `plotNodeProps()` draws the node-properties plane: node purity against
  cell count, joined with the ambiguity columns when provided. In the app it
  is a selectable plot box -- a lasso or click there selects nodes everywhere
  -- and when labels are present it offers selection by purity threshold.
- With labels present, the SOM 2D plots offer purity and the node's majority
  label among their color and size choices.
- `vignettes/spectral_flow_case_detailed_study.Rmd` gains a closing section
  running these functions on the case study's anchor SOM: the known gating
  plugs in as one colData column, and the section shows the bridge from the
  document's global pipeline metrics to the per-node interactive view.

## Bug fixes

* `R CMD check`'s "R code for possible problems" note is back to its
  pre-feature set: the names that `.applyLabelStatsToBundle()` restores from
  the saved bundle (`dList`, `dend_obj`, `dend_table`) are declared in
  `utils::globalVariables()`, like the other bundle names before them.

# CySA 0.99.40

## New features

- The ridge plot now sits directly above the SOM 2D plots, and clicking an
  outlier dot in the ridge plot steers the main 2D plot to that dot's
  marker: when the marker is not already shown, the first preset pair that
  contains it is activated (keeping the pair's own orientation), otherwise
  the marker is placed on the x axis with the first other marker from the
  marker list on the y axis. The clicked cell is then highlighted in the
  main 2D plot as a black diamond at its raw values on the active axes,
  together with a ring around its SOM node and a dashed connector between
  the two (in the plain, non-projection view) -- so outlier cells identified
  in the ridge plot can be located in the node-level 2D plot without
  manually cross-referencing axes. Click payloads that carry no marker
  (e.g. the bare-cell overlay trace, or an empty event) are ignored without
  error, and the highlight degrades gracefully when a marker has no SOM
  code column.

# CySA 0.99.39

## New features

- Added `vignettes/spectral_flow_case_study.Rmd`, a self-contained, reproducible
  case study demonstrating CySA on real spectral flow cytometry data: 56 FCS
  files (8 manually-gated populations x 7 donors, 24-marker Cytek Aurora panel)
  from an openly-downloadable Zenodo record (DOI 10.5281/zenodo.15723074), with
  the manual-gating labels hidden from FlowSOM clustering and used afterwards
  only to score the result (adjusted Rand index) and cross-check the new
  cluster-ambiguity diagnostic. Like `catalyst_workflow.Rmd`, it depends on an
  external download and is excluded from the Bioconductor-submission vignette
  build (`.Rbuildignore`). Added `flowCore` and `mclust` to `Suggests`.
- `spectral_flow_case_study.Rmd` now also runs `CyEt::optimizeBvalueMatrix()`
  (an independent, Shapiro-Wilk-normality-based per-marker arcsinh cofactor
  estimator) over all 24 markers and compares its choices against the
  vignette's own paper-derived cofactors, with cached results
  (`vignettes/support/cyet/`) so the ~2-hour, 24-marker optimization is not
  re-run on every vignette build. The comparison answers three questions
  directly: (1) the vignette's cofactors are already close to each marker's
  own ARI-optimal sweep value, not just a plausible guess; (2) CyEt's
  normality-optimal cofactors do *not* improve clustering here (ARI
  0.87-0.88 vs. 0.97-0.98) because normality-optimization pulls a bimodal
  marker's two populations together, which is the opposite of what
  SOM-based clustering wants; and (3) CyEt's own fits show a real, roughly
  two-fold systematic difference between type/lineage markers (median 368,
  excluding CD45) and state/functional markers (median 685) -- though this
  asymmetry has no effect on `cluster()`'s output, since `CATALYST::cluster()`
  only ever uses `type_markers(sce)` as SOM features. Also documents a
  genuine instability found in CyEt's own search for CD45 specifically (two
  different local optima depending on the tested cofactor range). Added
  `CyEt` to `Suggests`.
- `spectral_flow_case_study.Rmd` now also runs CyEt at the granularity it
  is meant to be used at -- per donor *and* per marker (168
  `optimizeBvalueVector()` calls, 6,000-cell-per-donor subsample, cached in
  `vignettes/support/cyet/cyet_per_donor_resolved.rds`) -- rather than
  pooling all donors into one subsample as the previous entry's comparison
  did. This surfaces real donor-to-donor cofactor heterogeneity the pooled
  run averaged away, and a second boundary artifact (GZMB, one donor) the
  pooled run also hid.
- `spectral_flow_case_study.Rmd` now tests the full three-step pipeline
  this enables: transform each donor with its own per-donor/per-marker
  CyEt cofactor, realign the donors with `cyCombine::batch_correct()`
  (per-cell-state-group `ComBat`), then an explicit per-marker z-score,
  all before FlowSOM, with cached results
  (`vignettes/support/cycombine/cycombine_per_donor_compare.rds`). None of
  the three steps helps: per-donor CyEt transformation alone drops ARI to
  0.634 (worse than the pooled-CyEt result), cyCombine's realignment
  recovers most but not all of that (0.823), and the z-score on top barely
  changes it further (0.821) -- all well short of this vignette's simple
  baseline of one shared per-marker cofactor with no correction (0.968).
  `cyCombine` is a GitHub-only package (not on CRAN/Bioconductor), so
  unlike `CyEt` it is deliberately **not** added to `Suggests`.
- `spectral_flow_case_study.Rmd` now documents, with a reproducible
  comparison (`vignettes/support/metaclustering/metaclustering_compare.rds`),
  that `CATALYST::cluster()`'s SOM-building step (`FlowSOM::ReadInput()` +
  `FlowSOM::BuildSOM()`) is genuinely FlowSOM's own code, while its
  metaclustering step uses `ConsensusClusterPlus::ConsensusClusterPlus()`
  rather than FlowSOM's own `metaClustering_consensus()` -- re-running
  metaclustering with FlowSOM's own function on the identical SOM codes
  gives ARI 0.977 vs. 0.968 for CATALYST's default, a real but small
  difference between two distinct consensus-clustering implementations.
- `spectral_flow_case_detailed_study.Rmd`'s `cyet_donor` cofactor set is
  regenerated with `CyEt` 0.99.13 (up from 0.99.12), which added adaptive
  range extension and explicit convergence status to
  `optimizeBvalueVector()`; this fixed a boundary artifact on GZMB (one
  donor) but left CD45 non-convergent on 5 of 7 donors, now resolved with
  the new version's bimodality-aware `detect_bimodal_split()` instead of
  an unconditional replace-with-global-median. CD45 is also now excluded
  from this vignette's SOM clustering features (`type_m`) -- while
  remaining in the cofactor/transform/batch-correction path -- because
  these FCS files are already manually gated into CD45+ leukocyte
  populations upstream, so CD45 has no negative population left to
  discriminate (confirmed directly: raw-scale minimum 8,820 across
  725,748 cells, zero CD45-negative cells) and its cofactor is for that
  reason only weakly identifiable. Re-running the full factorial after
  both changes shows per-donor CyEt remains, by a wide and stable margin,
  the one cofactor strategy clearly worse than the paper-derived cofactors
  (paired dARI -0.080, vs. -0.091 before either change) -- neither change
  alters that conclusion. Excluding CD45 from clustering is a small net
  cost in mean ARI across every cofactor strategy (`paper` 0.824 -> 0.803,
  best flat 0.851 -> 0.845), not an improvement, and is kept for the
  conceptual reason above rather than because it helps the numbers.

## Bug fixes

- `spectral_flow_case_study.Rmd` compared SOM cluster assignments against
  the manual-gating ground truth using the raw factor levels of
  `cluster_ids(sce, "som100")` instead of the underlying integer node IDs,
  which silently desynchronized the two whenever a factor's levels were not
  already sorted numerically; all ARI and cluster-ambiguity cross-checks in
  the vignette now go through
  `as.integer(as.character(cluster_ids(sce, "som100")))` consistently.
- Corrected several per-marker arcsinh cofactors in the vignette's panel
  definition that did not match the values implied by the source paper's own
  figures, and added a cofactor-sensitivity sweep (ARI vs. a range of flat
  cofactors) to justify the chosen values empirically rather than by
  inspection alone.

# CySA 0.99.38

## New features

* Added a "Cluster ambiguity" tab to the Sample diagnostics panel: for
  every cell, compares the distance to its own SOM node's code vector
  against the distance to the nearest *other* node, and bar-plots the
  percentage of each cluster's cells where those two distances are
  nearly tied (`computeClusterAmbiguity()`/`plotClusterAmbiguity()`).
  A high bar flags a node that is not well separated from its
  neighbours -- a scale-free, quantitative symptom of clustering on too
  many or poorly-discriminating marker dimensions, confirmed empirically
  on the Bodenmiller BCR-XL dataset (mean cell-to-code distance nearly
  doubled, and the fraction of ambiguous cells rose from 30% to 49%,
  when inflating an 8-marker lineage-only clustering to 17 markers).
  Restricted to whichever clusters/groups are selected under "Clusters/
  groups to include" -- picking a subset narrows the rival nodes
  considered for ambiguity, not just the cells shown.

# CySA 0.99.37

## New features

* Added a browser-tab favicon: four blue SOM nodes around one red
  "selected" node, the same red/blue convention already used throughout
  the app (e.g. the interactive dendrogram and 2D plots). The icon is an
  inline SVG data URI, so it needs no extra package resource file or
  `shiny::addResourcePath()` registration.

# CySA 0.99.36

## New features

* Added a "Sample diagnostics" panel reproducing the sample-level QC
  figures of the CATALYST/CyTOF reference workflow
  (`vignettes/catalyst_workflow.Rmd`) inside CySA itself: a pseudobulk MDS
  plot of the samples (equivalent of `CATALYST::pbMDS()`), a per-sample
  marker-expression heatmap (`plotExprHeatmap()`), and marker
  non-redundancy scores (`plotNRS()`), each optionally coloured by an
  `experiment_info` column.
* Added an "Ordered cluster heatmap" panel: marker medians per SOM node,
  ordered by the existing dendrogram's leaf order, with a cell-count bar
  alongside -- CySA's equivalent of `CATALYST::plotMultiHeatmap(bars =
  TRUE)`, reusing the dendrogram already built for the interactive
  dendrogram panel instead of a second, independently ordered clustering.
* Added a "Differential testing (diffcyt)" panel: GLMM-based differential
  abundance (DA) and LMM-based differential state (DS) tests, adapted
  from diffcyt's low-level pipeline (`calcCounts()`/`calcMedians()` +
  `testDA_GLMM()`/`testDS_LMM()`). Unlike the reference workflow, which
  tests CATALYST metaclusters, this tests CySA's own saved node groups
  (`outputList`), so the comparison of interest is always the grouping
  the user has already curated interactively. The UI enforces the two
  known diffcyt/GLMM gotchas found while building
  `vignettes/catalyst_workflow.Rmd`: a DA test requires a random effect
  (normally `sample_id`), and a DS test must not use `sample_id` as its
  random effect (each cluster-sample combination already contributes one
  median value, which perfectly confounds a `sample_id` random intercept
  with the residual). Adds `diffcyt` as a new package dependency.
* Every dashboard box now has a small info button next to its title that
  opens a modal explaining what the box does, what to expect, how to use
  it, and what it requires -- longer-form help than the one-line hover
  tooltips already on individual controls.
* "Sample diagnostics", "Ordered cluster heatmap" and "Differential
  testing" now let you restrict the analysis to a chosen selection of
  saved clusters/metaclusters (and, for the ordered heatmap, markers),
  instead of always using every cell/cluster. The differential-testing
  panel's own group selection replaces always testing every saved group.

## Bug fixes

* The differential-testing panel crashed with "NA/NaN/Inf in foreign
  function call" whenever the compared groups shared a SOM node -- most
  commonly, leaving the always-present "selected" pseudo-group selected
  alongside the named group just created from it. `.nodeGroupMap()`'s
  node-to-group lookup silently kept only the first group for a shared
  node, leaving the other group with zero cells; that produced an
  all-`NA` row in `diffcyt::calcMediansByClusterMarker()`'s output, which
  crashed `diffcyt::plotHeatmap()`'s internal `hclust()` row clustering.
  The compared groups are now checked up front and rejected with a clear
  error if they overlap, and a genuinely empty saved group (e.g. an
  emptied-out "Rest") is now dropped with a warning instead of reaching
  diffcyt at all. The default group selection also now excludes
  "selected", the most common source of this overlap.

# CySA 0.99.35

## User-interface improvements

* New optional neighbourhood-agreement view for the main SOM 2D plot
  (`somDataMain`): a "Neighbourhood agreement" checkbox (default off)
  colours each SOM node by how well the displayed marker pair preserves
  its nearest neighbours in the full SOM code space (Jaccard index of
  the node's k nearest neighbours in both spaces, k = 10). Nodes
  rendered blue are faithful to the full map; hot colours flag nodes
  whose local neighbourhood the pair distorts. The tooltip gains an
  "agreement:" line. Saved-group colouring still wins when active, and
  the checkbox can be combined with "Colour groups". Only SOM nodes
  with complete pair coordinates are scored.

## Bug fixes

* Naming the current selection ("Apply name") saved whatever node ids
  were typed into "Cluster numbers" verbatim, without checking them
  against the actual SOM. A stale or mistyped id (e.g. one left over
  after switching datasets, or a typo) saved into a named group later
  crashed every output that read that group by indexing
  `clusterPatientTable`/the SOM codes with it and hitting "subscript
  out of bounds": the counts and percentage bar plots, the cell-count
  table, the per-sample percentages, and both violin plots. Node ids
  are now validated against the current data when a name is applied,
  and every remaining place that indexes by a saved group's node ids
  (`.computeRelativeCounts()`, `.cellCountsTable()`, `.countComparison()`,
  `.violinData()`) additionally drops any id outside the current range
  instead of crashing, so a stale group degrades gracefully.
* `cellPercentages` threw "subscript out of bounds" whenever "Parent
  population for percentages" was set to a numeric `experiment_info`
  column (e.g. a cytometer's total-cell count): the sample-id row names
  needed to look up that column were built inline and lost by
  `apply()`/`as.numeric()`, instead of reusing the shared
  `.numericExperimentInfo()` helper that re-attaches them.

# CySA 0.99.34

## User-interface improvements

* The marker ridge plot is rebuilt on top of `ggridges`: a new geom,
  `geom_ridgeline_traced()`, subclasses `ggridges::GeomRidgeline()` and
  understands one additional value of the `datatype` aesthetic, `"path"`,
  so the ridges, the connected trace polylines and the outlier dots are
  drawn inside one layer. The traces therefore sit correctly above every
  ridge regardless of draw order, and the layer count drops from four to
  one.
* The marker ridge box in the app is now interactive (plotly,
  `plotRidgeInteractive()`). Each marker row shows the all-cell density,
  the density of the selected cells, a shaded band holding the central
  fraction of the selected cells ("Band coverage", default 0.99: the
  0.5% to 99.5% quantiles, so 1% of the selected cells fall outside on
  each marker), the selection mean (solid) and the all-cell mean
  (dotted). Selected cells outside the band are outlier dots; hovering
  shows the cell, node, sample, value and within-selection z-score.
  Clicking a dot traces that cell across every row and a new panel
  explains its SOM node membership: the per-marker contributions to the
  squared-distance difference between its own node and the nearest node
  outside the selection (clustering markers from `map$colsUsed`).
  Rows are sorted by separation from all cells (checkbox), clustering
  markers are marked with `*`, and markers on very different scales get a
  per-row scale. The static ggplot version stays in the PDF export; it
  converts the band coverage to the matching sd multiple (normal
  quantile, 0.99 -> 2.58 sd).

## Documentation

* The vignette and README gained a "Marker ridge: outliers and node
  membership" section with a screenshot
  (`man/figures/README-marker-ridge.png`) explaining the band, the
  outlier dots and the node membership panel.
* The vignette's table of statistical features, previously flattened
  into one paragraph, is a proper Markdown table again.

## Bug fixes

* The marker ridge plot collapsed into a single flat line at the top of
  the panel: `GeomRidgelineTraced` gave trace rows a draw-order `ymin`
  of a fixed 1 data unit below the ridges, which trained the y scale and
  dwarfed the (small) density values. The offset is now a tiny fraction
  of the data span.
* Marker row labels in the ridge plot were drawn outside the panel and
  clipped; they now sit just inside the left edge of each row.
* The static (PDF) marker ridge plot drew traced cells at the values of
  other cells whenever not every selected cell was extreme: trace values
  were looked up by the cell's position among the extreme cells instead of
  among all selected cells.
* The interactive ridge plot's outlier dots showed another marker's value
  when rows were sorted by separation (values were read by row position
  instead of marker name).
* Selection observers registered in a loop (SOM main plot, t-SNE, UMAP,
  PCA) now force their source id, so each listens to its own plot.
* Selecting many SOM nodes crashed the app ("C stack usage too close to
  the limit"): the ridge cluster picker pre-selected every selected node,
  drawing one density line per node and row (thousands of plotly traces).
  The picker now pre-selects the five largest nodes and accepts at most 20,
  and both ridge plots cap individual node lines at 20.
* The ridge cluster picker offered the clusters of the cells whose column
  index equalled a selected node id instead of the selected nodes.
* The "Rest" group vanished after the first group was named or removed:
  `.rebuildOutputList()` counted the previous Rest as used, so the new Rest
  was always empty.
* The FlowSOM marker pie output was sized for all selected nodes (up to
  8000 px high) although at most "Max pies" pies are drawn.

## Internal changes

* The Shiny server is restructured: all per-session state (static data,
  reactive values, reactives) lives in one environment built in stages by
  small builders (`.initServerState()`, `.build*Reactives()`), and the
  observer and output registrars are split into per-feature functions
  that take this state instead of ~50 positional arguments. The state is
  still exposed by name inside the server function, so `testServer()`
  expressions keep working.
* Every function is now at most 50 lines (BiocCheck recommendation), all R
  sources and vignette code use 4-space indents, and no line exceeds 80
  characters. Long plot, UI and preparation functions were split into
  documented helpers; repeated code (violin data, tooltip wrappers, zoom
  limits, customdata selection observers) is shared.
* Blanket `suppressWarnings()` / `suppressMessages()` calls are replaced by
  targeted handling: `.asIntegerQuiet()` parses node ids without coercion
  warnings, and the plotly conversion muffles only its known harmless
  warnings and its messages.
* PDF export draws plots with `methods::show()` instead of `print()`;
  `prepClusterSelector()` builds the dendrogram table without `<<-` and has
  a runnable example; unused pre-computations (per-channel limits, raster
  plot) were removed from the app builder.

## Bioconductor checks

* `R CMD check`: 0 errors, 0 warnings, 0 notes. Fixed: `ggridges` is now
  imported explicitly (`importFrom(ggridges, GeomRidgeline)`), the aes
  column names of `plotRidgeTraceFunc()` are declared in
  `utils::globalVariables()`, the duplicated `GeomRidgelineTraced` Rd
  alias is merged into `geom_ridgeline_traced.Rd`, `intro.html` plus the
  non-vignette files in `vignettes/` are excluded via `.Rbuildignore`,
  and the two remaining literal non-ASCII characters in `R/` are now
  `\uXXXX` escapes.
* `BiocCheck::BiocCheck()`: 0 errors, 0 warnings. The only remaining NOTEs
  are the roxygen2/`R CMD build`-generated `NAMESPACE` and `\usage{}`
  indentation (not hand-editable), the optional `fnd` Authors@R role, an
  optional biocViews suggestion, and the Bioc-Devel mailing list check
  (requires admin credentials to verify).

# CySA 0.99.33

## User-interface improvements

* The marker ridge plot no longer draws the connecting trace lines between
  marker rows: the outlier dots (one per traced cell and marker, sitting on
  the cluster ridge) remain, but the polylines cluttered the joyplot.
  `plotRidgeTraceFunc()` gained a `showTraces` argument (default `TRUE`)
  to bring the lines back if wanted.

# CySA 0.99.32

## Bug fixes

* Resizing plots no longer loops endlessly. `shinyjqui::jqui_resizable()`
  wraps outputs sized 100% of their parent; for outputs hidden at wrap
  time (collapsed boxes such as "marker ridge") the parent height is
  auto, making the wrapper content-driven: each render re-reported a new
  size and re-rendered the plot. Plot outputs are now pinned to a fixed
  pixel height by CSS, the window-resize script no longer re-imposes
  100% heights (width only), and the ridge plot is rendered at a fixed
  resolution so size changes are handled client-side instead of
  triggering server re-renders.

# CySA 0.99.31

## User-interface improvements

* The marker ridge plot now renders per-cluster densities as stacked
  joyplot ridges: each selected cluster's filled density curve is offset
  above the row base (classic ridgeline style), with the thicker black
  "combined" line sitting on the base of the row. Row heights account for
  the cluster-stack height so bands never overlap the row above.

## Bug fixes

* Trace vertices in the marker ridge plot were looked up by exact match on
  the density grid; a cell's value almost never falls exactly on a grid
  point, so every trace collapsed flat at the row base. Trace vertices are
  now interpolated (`stats::approx`) at the cell's value and sit on their
  cluster's density curve.
* The version tag in the dashboard header ("CySA x.y.z") is now rendered on
  its own line: it was previously appended inline after the title and got
  clipped off by the fixed-width AdminLTE logo box (`overflow: hidden`).
  The logo box also relaxes its inherited 50px line-height so both the
  title and the version line fit inside the 50px header.
* Cluster density lines and traces in the marker ridge plot now use whole
  Set2 brewer colours instead of interpolated ones: `colorRampPalette`
  produced washed-out greys for small numbers of clusters (e.g. grey
  #B3B3B3 for the second cluster), making clusters and traces nearly
  indistinguishable from the background.

# CySA 0.99.30

## User-interface improvements

* The marker ridge plot now distinguishes selected clusters: each selected
  cluster gets its own solid density line (capped at five; with more, the
  user picks which clusters get individual lines through a new picker that
  appears only in that case), plus a thicker black "combined" line pooling
  all selected clusters. The pastel background remains the full-dataset
  density.
* Trace lines and outlier dots are coloured by the traced cell's cluster,
  so multi-cluster selections remain distinguishable; the combined line is
  drawn thicker (1.1 vs 0.7) so the aggregate reads as one.
* Fixed duplicate ggplot colour scales that flooded the log with
  "Adding another scale for colour" warnings and silently dropped
  selected-cluster density rows; all ridge layers now share one unified
  colour scale.

# CySA 0.99.29

## User-interface improvements

* The dashboard header now shows the running package version
  (e.g. "Cluster Selector CySA 0.99.29") so users and bug reports always
  reference the right build.
* The marker ridge plot is debounced (500 ms) so dragging a resizable
  corner does not continuously re-render the plot while resizing.

# CySA 0.99.8

## User-interface improvements

* Added hover-over tooltips (`title` attributes) to every non-obvious control
  in the `clusterSelector()` sidebar and body panels.
* Added inline `helpText()` explanations for subtle controls such as selection
  mode, scatter auto-zoom percentile, and dimension-reduction settings.
* Added selection hints to all plot outputs that support click/lasso
  selection.
* Extended roxygen comments for all UI builder helpers in `R/ui-helpers.R`.
* Added a "marker ridge" box: one row per selected marker showing the density
  of all cells (filled pastel) and of the currently selected clusters (solid
  line), with dashed threshold lines at `mean +/- threshold * sd` and capped
  trace lines connecting extreme cells across the marker rows. Each traced
  cell is also marked with a dot at its value in every row, and traces are
  coloured by the marker on which a cell deviates most, making it easy to
  see why an outlying cell belongs to the selection. Rows carry headroom so
  the ridgeline bands stay separated.
* Fixed a crash in the marker ridge trace extraction when the assay matrix
  has missing or misaligned cell colnames: trace values are now read by
  position instead of by cell name.

# CySA 0.99.7

## Bioconductor submission fixes

* Converted vignette R code chunks to 4-space indentation.

# CySA 0.99.6

## Bioconductor submission fixes

* Ensured all R source indentation is an exact multiple of 4 spaces.

# CySA 0.99.5

## Bioconductor submission fixes

* Fixed remaining indentation outliers (hanging indents and manually edited
  blocks) so all R source lines use multiples of 4 spaces.

# CySA 0.99.4

## Bioconductor submission fixes

* Converted R source files from 2-space to 4-space indentation to conform to
  Bioconductor style guidelines.

# CySA 0.99.3

## Bioconductor submission fixes

* Wrapped the second occurrence of a long `outputList[["Rest"]]` assignment.

# CySA 0.99.2

## Bioconductor submission fixes

* Wrapped remaining long lines flagged by `BiocCheck`.

# CySA 0.99.1

## Bioconductor submission fixes

* Replaced `<<-` with `shiny::reactiveVal()` for app state.
* Removed debug `cat()`/`print()` calls; kept `print()` only inside PDF device
  context.
* Applied `styler` formatting and wrapped selected long lines.

# CySA 0.99.0

## New features

* Initial Bioconductor submission of `CySA`, an interactive Shiny dashboard for
  selecting and visualizing clusters from flow-cytometry data stored in
  `SingleCellExperiment` objects.
* `clusterSelector()` returns a Shiny app object that can be launched with
  `shiny::runApp()`.
* `prepClusterSelectorData()` subsamples a `SingleCellExperiment` and builds the
  inputs required by `clusterSelector()`.
* `plotSOMScatter()` and `plotCytoScatter()` provide ggplot2-based SOM and scatter
  visualizations.

## Bioconductor readiness

* Replaced broad package-level `@import` directives with targeted
  `@importFrom` roxygen tags to reduce namespace masking.
* Removed the archived `ggthemr` dependency; theming now uses `theme_minimal()`
  and explicit `scale_*_manual()` scales.
* Added `NEWS.md`, expanded `README.Rmd`, and a package vignette.
