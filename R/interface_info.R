# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# interface_info.R ----
# A small info ("i") button on every dashboard box, opening a modal that
# explains what the box does, what to expect, how to use it, and what it
# requires -- longer-form help than the one-line hover tooltips in
# .clusterSelectorTooltips (interface_clusterSelector.R). Each box passes
# its key into .boxInfoRegistry as .primaryBox(..., infoId = "<key>");
# .registerBoxInfoObservers() wires up every key in the registry once.

#' Help Text for Every Dashboard Box
#'
#' Named by the \code{infoId} passed to \code{.primaryBox()}. Each entry's
#' \code{title} becomes the modal's title and \code{body} a character
#' vector of paragraphs (rendered one \code{<p>} per element).
#'
#' @keywords internal
.boxInfoRegistry <- list(
    interactiveTree = list(
        title = "Interactive tree",
        body = c(
            paste(
                "What it does: shows the saved cluster groups as a",
                "collapsible tree, so you can see and navigate the",
                "current grouping hierarchy at a glance."
            ),
            paste(
                "What to expect: only appears when the app was launched",
                "with a colTree object; otherwise this box is absent."
            ),
            paste(
                "How to use it: expand/collapse branches to explore",
                "groups; this is a navigation aid, it does not itself",
                "change the selection."
            ),
            paste(
                "Requires: nothing from you -- it reflects whatever",
                "groups already exist in the sidebar's \"Named groups\"",
                "list."
            )
        )
    ),
    twoDPlot = list(
        title = "2D plot",
        body = c(
            paste(
                "What it does: the main interactive scatter plot of all",
                "SOM nodes on two chosen marker axes (set in the",
                "\"SOM 2D plots\" box below)."
            ),
            paste(
                "What to expect: one point per SOM node, positioned by",
                "its marker expression; colour/size reflect whatever",
                "statistic you picked there."
            ),
            paste(
                "How to use it: click a point, or drag a lasso/box, to",
                "add those nodes to the current selection. The",
                "sidebar's \"Selection mode\" (view/add/remove/remove",
                "others) controls how a new click combines with the",
                "existing selection."
            ),
            paste(
                "Requires: nothing extra -- it is always available",
                "once the app has loaded."
            )
        )
    ),
    dendrogram = list(
        title = "Interactive dendrogram",
        body = c(
            paste(
                "What it does: a hierarchical clustering of the SOM",
                "nodes (by marker-profile similarity), drawn as a",
                "dendrogram you can click on."
            ),
            paste(
                "What to expect: leaves are individual SOM nodes;",
                "selected nodes are highlighted in red among the",
                "others in blue."
            ),
            paste(
                "How to use it: click a leaf (single node) or an",
                "internal branch (every node under it) to select those",
                "nodes, following the same \"Selection mode\" rule as",
                "the 2D plot."
            ),
            paste(
                "Requires: a precomputed dendrogram (the",
                "\"dend\"/\"dendTable\" arguments passed to",
                "clusterSelector())."
            )
        )
    ),
    som2dPlots = list(
        title = "SOM 2D plots",
        body = c(
            paste(
                "What it does: lets you pick which two markers (or a",
                "preset pair) define the axes of the main 2D plot",
                "above, plus optional colouring, and shows",
                "t-SNE/UMAP/PCA projections of the SOM nodes for a",
                "complementary view."
            ),
            paste(
                "What to expect: changing the X/Y axis selectors",
                "immediately redraws the main 2D plot; the three small",
                "plots at the bottom (t-SNE, UMAP, PCA) are computed",
                "from the SOM code vectors using whichever markers you",
                "pick under \"Markers for dim. reduction\", not the",
                "X/Y axes."
            ),
            paste(
                "How to use it: choose axes directly, or pick a saved",
                "preset pair; toggle \"Colour groups\" to colour by",
                "saved group instead of by statistic, or",
                "\"Neighbourhood agreement\" to flag nodes whose",
                "neighbours this marker pair distorts. All of the",
                "small plots support the same click/lasso selection as",
                "the main 2D plot."
            ),
            paste(
                "Requires: nothing extra for the main plot;",
                "t-SNE/UMAP/PCA need at least one marker selected",
                "under \"Markers for dim. reduction\" and recompute",
                "whenever that selection, the perplexity, or",
                "n_neighbors changes."
            )
        )
    ),
    staticSomViews = list(
        title = "All SOM pair views (static)",
        body = c(
            paste(
                "What it does: a grid of extra, non-interactive SOM",
                "scatter plots, one per marker pair you choose to keep",
                "around -- a quick way to scan many marker",
                "combinations at once without repeatedly changing the",
                "main plot's axes."
            ),
            paste(
                "What to expect: each small plot marks the currently",
                "selected SOM nodes as red dots, but (unlike the main",
                "2D plot) you cannot click on these to change the",
                "selection."
            ),
            paste(
                "How to use it: pick an X/Y pair and click \"Add\" to",
                "append it to the grid; pick an existing pair from",
                "\"Remove pair\" and click the trash icon to drop it.",
                "\"Columns per row\" and \"Plot height\" control the",
                "grid layout."
            ),
            paste(
                "Requires: expand the box to actually render the plots",
                "(they are not drawn while collapsed, to save time",
                "when you don't need them)."
            )
        )
    ),
    flowSomPies = list(
        title = "FlowSOM marker pies",
        body = c(
            paste(
                "What it does: draws one pie chart per currently",
                "selected SOM node, with each slice sized by that",
                "node's mean expression of one marker -- a compact",
                "per-node marker-profile summary."
            ),
            paste(
                "What to expect: pies use a consistent",
                "colour-to-marker mapping across nodes, so the same",
                "slice colour always means the same marker; only the",
                "first \"Max pies\" selected nodes are drawn."
            ),
            paste(
                "How to use it: adjust \"Columns\", \"Size per pie\"",
                "and \"Max pies\" to control the layout; select/",
                "deselect SOM nodes elsewhere in the app (2D plot,",
                "dendrogram, raster, ...) to change which pies are",
                "shown."
            ),
            paste(
                "Requires: at least one SOM node selected; nothing is",
                "drawn otherwise."
            )
        )
    ),
    flowSomStars = list(
        title = "FlowSOM stars",
        body = c(
            paste(
                "What it does: an interactive FlowSOM \"star\" plot --",
                "each node is a star whose spokes show marker",
                "expression, laid out on FlowSOM's own",
                "minimum-spanning-tree grid and coloured by",
                "metaclustering."
            ),
            paste(
                "What to expect: this box only appears when the app",
                "was launched with a FlowSOM object (the \"fsom\"",
                "argument); otherwise it is absent entirely."
            ),
            paste(
                "How to use it: click or lasso nodes to add them to",
                "the current selection, exactly like the main 2D plot."
            ),
            paste(
                "Requires: a FlowSOM object passed to clusterSelector();",
                "nothing further from you."
            )
        )
    ),
    stats = list(
        title = "Stats",
        body = c(
            paste(
                "What it does: tabulates and plots the cell counts of",
                "the currently selected SOM nodes, as raw counts, as a",
                "percentage of a reference population, and (below, in",
                "\"t-Test\") compared between two groups of samples."
            ),
            paste(
                "What to expect: the table/bar plots update as soon as",
                "you change the selection elsewhere; \"Parent",
                "population for percentages\" changes only the",
                "denominator used for the percentage column and bar",
                "plot, not the counts table."
            ),
            paste(
                "How to use it: pick a numeric experiment_info column",
                "or saved group under \"select comparison Stats\" to",
                "add a side-by-side comparison series to the count bar",
                "plot; pick a parent population (another numeric",
                "column, or \"none\" for the grand total) for the",
                "percentage bar plot; \"single SOM node counts\" looks",
                "up one specific node's counts regardless of the",
                "current selection."
            ),
            paste(
                "Requires: at least one SOM node selected for the",
                "counts/percentages to be meaningful."
            )
        )
    ),
    tTest = list(
        title = "t-Test",
        body = c(
            paste(
                "What it does: runs a two-sample t-test comparing the",
                "(optionally percentage-normalised) cell counts of the",
                "current selection between two groups of samples that",
                "you define."
            ),
            paste(
                "What to expect: a plain-text result (t statistic,",
                "degrees of freedom, p-value); it prints",
                "\"not enough data\" when both groups have zero",
                "variance (e.g. too few samples)."
            ),
            paste(
                "How to use it: pick a factor experiment_info column",
                "under \"Select groups variable for t-test\" (this",
                "populates the group-1/group-2 choices), then assign",
                "samples to \"group 1\" and \"group 2\". The same",
                "\"Parent population for percentages\" from the Stats",
                "box above applies here too."
            ),
            paste(
                "Requires: at least one SOM node selected, a factor",
                "experiment_info column, and at least 2 samples per",
                "group with non-zero variance in their counts."
            )
        )
    ),
    somRaster = list(
        title = "som plots",
        body = c(
            paste(
                "What it does: a faceted raster heatmap of every SOM",
                "node's marker expression on the SOM's own 2D grid",
                "layout (one facet per marker), with a companion plot",
                "for clicking/lassoing nodes directly on that same",
                "grid."
            ),
            paste(
                "What to expect: the heatmap facets are static",
                "(colour = expression, viridis scale); the lower",
                "plotly grid highlights currently selected nodes and",
                "lets you change the selection."
            ),
            paste(
                "How to use it: click or lasso nodes in the lower grid",
                "plot, following the sidebar's \"Selection mode\"; use",
                "the heatmap facets above it to see which markers a",
                "region of the grid expresses highly."
            ),
            paste(
                "Requires: a somRasterData argument (and optionally",
                "somRasterObj) passed to clusterSelector()."
            )
        )
    ),
    violinPlots = list(
        title = "violin plots",
        body = c(
            paste(
                "What it does: two violin plots comparing",
                "marker-expression distributions across your saved",
                "groups -- the first pools all SOM-code values per",
                "group, the second uses the raw per-cell data of the",
                "subsampled SCE."
            ),
            paste(
                "What to expect: groups with fewer than 2 SOM nodes",
                "are dropped from the first plot (there is no",
                "distribution to show); both plots redraw whenever you",
                "add/rename/remove a saved group or change which",
                "groups are compared (via the UpSet box's \"Select",
                "groups to show\", which both plots reuse)."
            ),
            paste(
                "How to use it: pick which markers to show under",
                "\"Select markers to show\"; create/rename/remove",
                "groups from the sidebar to change what gets compared."
            ),
            paste(
                "Requires: at least 2 saved groups with at least 2 SOM",
                "nodes each to get a meaningful comparison."
            )
        )
    ),
    ridgeBox = list(
        title = "marker ridge",
        body = c(
            paste(
                "What it does: a ridgeline plot of every chosen",
                "marker's expression density, overlaying \"all cells\"",
                "against the current selection, plus a per-marker band",
                "covering a chosen fraction of the selection's cells;",
                "cells outside that band are shown as individually",
                "clickable outlier dots. Clicking an outlier dot opens",
                "a second plot explaining which SOM node that specific",
                "cell belongs to."
            ),
            paste(
                "What to expect: rows can be sorted by how far the",
                "selection's mean is shifted from the all-cell mean",
                "(in all-cell standard deviations), to surface the",
                "markers most affected by your selection first; when",
                "more than five clusters are selected, a picker lets",
                "you choose which ones get their own density line (the",
                "combined line always covers the whole selection)."
            ),
            paste(
                "How to use it: choose markers under",
                "\"Select markers to show\"; adjust \"Band coverage\"",
                "to widen/narrow what counts as an outlier, and",
                "\"Max outlier cells\" to cap how many outlier dots are",
                "drawn (the most extreme are kept); click any outlier",
                "dot to see its SOM-node explanation below."
            ),
            paste(
                "Requires: at least one SOM node selected."
            )
        )
    ),
    upsetPlot = list(
        title = "UpSet plot",
        body = c(
            paste(
                "What it does: an UpSet plot showing how your saved",
                "groups overlap in SOM-node and cell-count terms --",
                "which nodes/cells are unique to one group versus",
                "shared across several."
            ),
            paste(
                "What to expect: the top bars count nodes/cells per",
                "overlap combination; the right-hand bars count",
                "nodes/cells per selected group on its own; the groups",
                "chosen here are also reused by both violin plots."
            ),
            paste(
                "How to use it: pick 2 or more saved groups under",
                "\"Select groups to show\" (more than 31 groups or",
                "more than 100 distinct overlap combinations are",
                "rejected, to keep the plot readable -- narrow the",
                "selection if that happens)."
            ),
            paste(
                "Requires: at least 2 saved groups."
            )
        )
    ),
    diagnostics = list(
        title = "Sample diagnostics",
        body = c(
            paste(
                "What it does: sample-level QC, reproducing CATALYST's",
                "reference CyTOF-workflow figures inside CySA -- a",
                "pseudobulk MDS plot of the samples, a per-sample",
                "marker-expression heatmap, and per-marker",
                "non-redundancy scores (NRS); see",
                "vignettes/catalyst_workflow.Rmd for the original",
                "figures these are modelled on."
            ),
            paste(
                "What to expect: all three tabs are computed from",
                "per-sample pseudobulk (median marker expression per",
                "sample) restricted to whichever clusters/groups you",
                "pick below -- this is a sample-centric view (\"what",
                "do my samples look like, within this cell",
                "population?\"), not a per-cluster comparison (use",
                "\"Ordered cluster heatmap\" below for that)."
            ),
            paste(
                "How to use it: pick which saved clusters/groups to",
                "include (leave every group selected to use the whole",
                "dataset); optionally colour samples by an",
                "experiment_info column (e.g. treatment/batch) to spot",
                "separation or outlier samples."
            ),
            paste(
                "Requires: at least 3 samples for the MDS tab (fewer",
                "markers/samples still work for the heatmap and NRS",
                "tabs, just with less to compare)."
            )
        )
    ),
    orderedHeatmap = list(
        title = "Ordered cluster heatmap",
        body = c(
            paste(
                "What it does: a marker-by-cluster heatmap of the",
                "chosen clusters/groups, ordered by the existing",
                "dendrogram's leaf order, with a cell-count bar",
                "alongside -- CySA's equivalent of",
                "CATALYST::plotMultiHeatmap(bars = TRUE), reusing the",
                "same dendrogram as the \"Interactive dendrogram\" box",
                "instead of a second, independently ordered",
                "clustering."
            ),
            paste(
                "What to expect: rows are clusters/SOM nodes (in",
                "dendrogram order), columns are markers; this is the",
                "box to use for directly comparing marker profiles",
                "across multiple clusters/metaclusters side by side."
            ),
            paste(
                "How to use it: pick which clusters/groups to show as",
                "rows, and which markers to show as columns."
            ),
            paste(
                "Requires: at least one cluster/group and one marker",
                "selected."
            )
        )
    ),
    diffTest = list(
        title = "Differential testing (diffcyt)",
        body = c(
            paste(
                "What it does: diffcyt GLMM-based differential",
                "abundance (DA) and LMM-based differential state (DS)",
                "tests, comparing your chosen saved",
                "clusters/metaclusters across samples -- adapted from",
                "the reference diffcyt workflow in",
                "vignettes/catalyst_workflow.Rmd, but testing CySA's",
                "own saved groups instead of CATALYST metaclusters."
            ),
            paste(
                "What to expect: DA tests whether a cluster's relative",
                "abundance differs between the two levels of your",
                "fixed-effect column; DS tests whether a marker's",
                "median expression within a cluster does. The compared",
                "groups must not overlap in SOM nodes -- a common",
                "mistake is leaving \"selected\" selected alongside the",
                "group you just named from it, since they then share",
                "every node; CySA rejects this with a clear error",
                "instead of crashing. Fitting can take a few seconds."
            ),
            paste(
                "How to use it: pick >= 2 non-overlapping saved groups",
                "to compare, the analysis type, a 2-level",
                "experiment_info factor as the fixed effect (the",
                "comparison of interest), and a random effect (DA",
                "requires one, normally sample_id; DS must NOT use",
                "sample_id -- pick \"none\" or a column like",
                "patient_id). For DS, also pick which markers to test.",
                "Click \"Run test\"."
            ),
            paste(
                "Requires: >= 2 non-overlapping saved groups, a 2-level",
                "factor column in experiment_info, and (for DA) a",
                "random-effect column."
            )
        )
    )
)

#' Info ("i") Button Opening a Box's Help Modal
#'
#' @param infoId Key into \code{.boxInfoRegistry}.
#'
#' @return An \code{actionButton} tag.
#'
#' @keywords internal
.infoButton <- function(infoId) {
    shiny::actionButton(
        inputId = paste0("boxInfo_", infoId), label = NULL,
        icon = shiny::icon("info-circle"),
        class = "btn-xs btn-info", style = "float:right; margin-left:8px;",
        title = "About this panel"
    )
}

#' Register the Info-Button Modal Observers
#'
#' One \code{observeEvent} per entry of \code{.boxInfoRegistry}, each
#' showing that entry's help text in a modal when its info button is
#' clicked.
#'
#' @param input,session Shiny objects.
#'
#' @return Called for its side effects.
#'
#' @keywords internal
.registerBoxInfoObservers <- function(input, session) {
    lapply(names(.boxInfoRegistry), function(infoId) {
        local({
            id <- infoId
            info <- .boxInfoRegistry[[id]]
            shiny::observeEvent(input[[paste0("boxInfo_", id)]], {
                shiny::showModal(shiny::modalDialog(
                    title = info$title,
                    shiny::tagList(lapply(info$body, shiny::tags$p)),
                    easyClose = TRUE, size = "l"
                ))
            }, ignoreInit = TRUE)
        })
    })
    invisible(NULL)
}
