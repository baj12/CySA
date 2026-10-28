## refine_engine.R — engine for the hierarchical SOM refinement experiment.
##
## Companion of vignettes/refinement_experiment.Rmd (design doc:
## docs/superpowers/specs/2026-10-28-hierarchical-refinement-experiment-design.md).
##
## Provenance: the run engine (cached, save_rds_atomic, show_plot, ari_score,
## node_purity, node_distances, ambiguity_metrics, meta_cluster, eta2,
## apply_scaling, apply_weights, prepare_X, assign_nodes, run_one, run_stage,
## get_base, transform_raw, get_cofactors, parse_name) is copied from
## vignettes/spectral_flow_case_detailed_study.Rmd so both documents score
## runs identically. Refinement-specific code (amb_cells, rival_of,
## gap_scores, refine_tree, flat_mean_run, expert_bound, stage drivers) is
## new and marked "refine:".
##
## Usage:
##   - SLURM stage drivers: set out_dir/n_cores, then
##     source("vignettes/refine_engine.R"); stage_xxx()
##   - vignette: sets params and out_dir BEFORE sourcing; dev_n > 0
##     subsamples the data for smoke renders (all heavy stages are cache
##     hits then; fresh dev-mode computations use the subsample).

## ---------------------------------------------------------------- config --

suppressPackageStartupMessages({
    library(FlowSOM)
    library(matrixStats)
    library(mclust)
})

if (!exists("out_dir")) out_dir <- "spectral_analysis_refine"
if (!exists("n_cores")) n_cores <- 4L
if (!exists("dev_n"))   dev_n <- 0L
for (d in c("", "figures", "runs")) {
    dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

## Caches: refinement runs live in out_dir (spectral_analysis_refine/), but
## the shared 00_inputs.rds (raw matrix, populations, donors, full cohort)
## is read from the detailed study's cache when present -- read-only, never
## rewritten by this engine.
refine_inputs_file <- file.path(out_dir, "00_inputs.rds")
shared_inputs_file <- "spectral_analysis/00_inputs.rds"

## ------------------------------------------------------------- helpers ----

#' Atomically write an RDS file.
#' @param obj object to save
#' @param path destination path
save_rds_atomic <- function(obj, path) {
    tmp <- paste0(path, ".tmp", Sys.getpid())
    saveRDS(obj, tmp)
    file.rename(tmp, path)
    invisible(path)
}

#' Compute-once helper.
#' @param file file name relative to `out_dir` (should end in .rds)
#' @param expr expression producing the object; evaluated lazily, only when
#'   the cache file is missing
#' @param csv if the result is a data.frame, also write a .csv next to it
#' @return the object (read from disk or freshly computed)
cached <- function(file, expr, csv = TRUE) {
    path <- file.path(out_dir, file)
    if (file.exists(path)) {
        message("[cache hit]  ", path)
        return(readRDS(path))
    }
    message("[computing]  ", path)
    res <- expr
    save_rds_atomic(res, path)
    if (csv && is.data.frame(res)) {
        utils::write.csv(res, sub("\\.rds$", ".csv", path), row.names = FALSE)
    }
    res
}

#' Print a ggplot and save it as PNG in `out_dir/figures`.
#' Re-rendered on every call (cheap; unlike stage caches the figure must
#' follow dev-mode subsampling).
show_plot <- function(p, name, width = 8, height = 5) {
    f <- file.path(out_dir, "figures", paste0(name, ".png"))
    ggsave(f, p, width = width, height = height, dpi = 150)
    print(p)
    invisible(p)
}

dist_labels <- c(`1` = "manhattan", `2` = "euclidean", `3` = "chebyshev", `4` = "cosine")

## ----------------------------------------------------------------- data ----

#' Parse "<population>_donor<N>.fcs" file names.
parse_name <- function(fn) {
    base  <- sub("\\.fcs$", "", basename(fn))
    donor <- as.integer(sub("_donor", "", regmatches(base, regexpr("_donor[0-9]+$", base))))
    pop   <- gsub(" ", "_", sub("_donor[0-9]+$", "", base))
    list(pop = pop, donor = donor)
}

common_markers <- c(
    "CCR4", "GZMB", "HELIOS", "GNLY", "CCR7", "CD4", "CD56", "CD3", "CD19",
    "CD95", "CD8", "CD103", "CD16", "CD14", "CD69", "CCR6", "HLA-DR", "CD25",
    "CD45RA", "PD1", "CD39", "FOXP3", "TIGIT", "CD45"
)
state_m <- c("CD69", "CD25", "CD39", "PD1", "TIGIT", "FOXP3", "GZMB", "CCR4",
             "CCR6", "HELIOS")
all_m  <- common_markers
## CD45 excluded from the SOM features (already CD45+ gated upstream), as in
## the detailed study. HELIOS is a state marker and never clustered.
type_m <- setdiff(common_markers, c(state_m, "CD45"))
## refine: the d5-7 CCR7 channel artifact (AF488 pedestal; discussion_02)
## excludes CCR7 from the clustering features. HELIOS needs no action
## (already a state marker). 12 markers remain.
type_m12 <- setdiff(type_m, "CCR7")

if (isTRUE(getOption("refine_synthetic", FALSE))) {
    ## Synthetic mode (tests only): the caller provides raw/pop/donor and has
    ## NOT loaded the cohort. Everything downstream (dev subsample, amb_idx,
    ## cofactors) works unchanged on the caller's data.
    N_full <- nrow(raw)
    message("[engine] synthetic mode: caller-provided raw/pop/donor")
} else {
    if (file.exists(shared_inputs_file) && !file.exists(refine_inputs_file)) {
        message("[inputs] copying shared cache ", shared_inputs_file)
        file.copy(shared_inputs_file, refine_inputs_file)
    }

    inputs <- cached("00_inputs.rds", local({
        fcs_files <- list.files("spectral_fcs", pattern = "\\.fcs$", full.names = TRUE)
        stopifnot(length(fcs_files) == 56)
        parts <- lapply(fcs_files, function(fp) {
            info <- parse_name(fp)
            ff <- flowCore::read.FCS(fp, transformation = FALSE, truncate_max_range = FALSE)
            pd <- pData(parameters(ff))
            idx <- match(common_markers, as.character(pd$desc))
            stopifnot(!anyNA(idx))
            m <- exprs(ff)[, idx, drop = FALSE]
            colnames(m) <- common_markers
            storage.mode(m) <- "double"
            list(m = m, pop = rep(info$pop, nrow(m)), donor = rep(info$donor, nrow(m)))
        })
        list(
            raw        = do.call(rbind, lapply(parts, `[[`, "m")),   # raw intensities, cells x 24
            population = factor(unlist(lapply(parts, `[[`, "pop"))),
            donor      = as.integer(unlist(lapply(parts, `[[`, "donor")))
        )
    }), csv = FALSE)

    raw <- inputs$raw; donor <- inputs$donor
    pop <- factor(stringi::stri_trans_general(
        stringi::stri_trans_nfc(enc2utf8(as.character(inputs$population))),
        "Latin-ASCII"))
    N_full <- nrow(raw)
}
if (dev_n > 0 && dev_n < N_full) {
    set.seed(1)
    keep <- sort(sample(N_full, dev_n))
    raw <- raw[keep, ]; pop <- droplevels(pop[keep]); donor <- donor[keep]
}
N <- nrow(raw)
donor_lab <- paste0("donor", donor)
pop_int <- as.integer(pop)
cat(N, "cells,", ncol(raw), "markers,", length(type_m), "type markers (",
    length(type_m12), "after CCR7 exclusion)\n")

## Fixed ambiguity-evaluation subset (comparable across all configurations;
## identical convention to the detailed study: seed 99, min(n_amb, N)).
## In dev mode it is intersected with the kept cells.
n_amb <- 50000L
set.seed(99)
amb_idx_full <- sort(sample(N_full, min(n_amb, N_full)))
if (dev_n > 0 && dev_n < N_full) amb_idx_full <- amb_idx_full[amb_idx_full %in% keep]
amb_idx <- match(amb_idx_full, if (dev_n > 0 && dev_n < N_full) keep else seq_len(N_full))
stopifnot(!anyNA(amb_idx))

## Cofactor set: hand-tuned per-marker cofactors from Tursi et al. 2025
## (spectral_flow_case_study.Rmd), 6000 elsewhere.
spectral_cofactors <- setNames(rep(6000, length(all_m)), all_m)
spectral_cofactors[c("CD14", "CD19", "CD3", "CD56", "CD45RA", "CD8", "CD4", "CCR7")] <-
    c(10000, 1000, 3000, 2000, 4000, 3000, 5000, 6000)

#' arcsinh-transform the raw matrix (cells x 24) with a cofactor set.
transform_raw <- function(cof) {
    asinh(sweep(raw, 2, cof[all_m], "/"))
}

#' Cells x type-marker matrix after cofactor transform.
#' Batch correction (cyCombine) is not used by the refinement experiment
#' (anchor config: batch = none).
get_base <- function(cf = "paper", markers = type_m) {
    stopifnot(cf == "paper")
    tr <- transform_raw(get_cofactors(cf))
    tr[, markers, drop = FALSE]
}

get_cofactors <- function(name) {
    stopifnot(name == "paper")
    spectral_cofactors
}

## ------------------------------------------------------------- metrics ----

#' Adjusted Rand index between two partitions (mclust::adjustedRandIndex).
ari_score <- function(a, b) mclust::adjustedRandIndex(a, b)

#' Node purity: fraction of cells that belong to the majority known population
#' of the cluster they were mapped to. Deliberately NOT adjusted for the
#' number of clusters; compare at matched cluster counts.
node_purity <- function(node, pop) {
    tab <- table(node, pop)
    sum(apply(tab, 1, max)) / length(node)
}

#' Purity of one cluster's cells: majority-population fraction. refine:
node_purity_one <- function(cl_idx) {
    if (length(cl_idx) == 0L) return(NA_real_)
    max(tabulate(pop_int[cl_idx])) / length(cl_idx)
}

#' Distance of every cell to every code vector, for FlowSOM's `distf`.
#'   1 = Manhattan  sum_p |x_p - c_p|
#'   2 = Euclidean  sqrt(sum_p (x_p - c_p)^2)
#'   3 = Chebyshev  max_p |x_p - c_p|
#'   4 = cosine     1 - <x, c> / (||x|| ||c||)
node_distances <- function(X, codes, distf) {
    n <- nrow(X); K <- nrow(codes)
    if (distf == 4) {
        eps <- .Machine$double.eps
        nx <- pmax(sqrt(rowSums(X^2)), eps)
        nc <- pmax(sqrt(rowSums(codes^2)), eps)
        return(1 - (X / nx) %*% t(codes / nc))
    }
    D <- matrix(0, n, K)
    for (j in seq_len(K)) {
        diff <- X - rep(codes[j, ], each = n)
        D[, j] <- switch(as.character(distf),
            "1" = rowSums(abs(diff)),
            "2" = sqrt(rowSums(diff^2)),
            "3" = rowMaxs(abs(diff)),
            stop("unsupported distf"))
    }
    D
}

#' Cluster ambiguity (CySA's "Cluster ambiguity" tab): own/rival distance
#' ratio per cell; ratio > 0.9 = ambiguous.
ambiguity_metrics <- function(X, codes, node, distf, idx) {
    D <- node_distances(X[idx, , drop = FALSE], codes, distf)
    own_j <- node[idx]
    n <- length(idx)
    own <- D[cbind(seq_len(n), own_j)]
    D[cbind(seq_len(n), own_j)] <- Inf
    rival <- rowMins(D)
    ratio <- own / pmax(rival, .Machine$double.eps)
    list(pct_ambiguous = 100 * mean(ratio > 0.9),
         median_ratio = median(ratio))
}

#' FlowSOM consensus metaclustering of code vectors (tempdir: no file
#' collisions between parallel workers).
#'
#' Upstream bug workaround: ConsensusClusterPlus's diagnostic tracking plot
#' (clusterTrackingPlot) always crashes at maxK = 2 (colorM has one row, so
#' matrix subsetting drops dims and 1:nrow(m) gets NULL) -- verified for code
#' matrices of 9/16/25/100 rows; k >= 3 is unaffected. The clustering itself
#' is computed identically for every k; only the diagnostic plot is skipped,
#' so the patch does not change any result. Patched once at engine load.
assignInNamespace("clusterTrackingPlot", function(m) invisible(NULL),
                   ns = "ConsensusClusterPlus")

meta_cluster <- function(codes, k, seed) {
    old <- setwd(tempdir()); on.exit(setwd(old), add = TRUE)
    FlowSOM::metaClustering_consensus(codes, k = k, seed = seed)
}

#' Between-population variance fraction (eta^2) of every column of X.
#' ORACLE: uses the known labels (diagnostic / upper bound only).
eta2 <- function(X, grp) {
    gs <- rowsum(X, grp)
    ng <- as.numeric(table(grp)[rownames(gs)])
    gm <- gs / ng
    mu <- colMeans(X)
    ss_between <- colSums(ng * sweep(gm, 2, mu, "-")^2)
    ss_total   <- colSums(X^2) - nrow(X) * mu^2
    ss_between / ss_total
}

## refine: chunked ambiguity for an ARBITRARY partition. Unlike
## ambiguity_metrics() this takes cluster *representatives* (rows of R, e.g.
## cell-mean profiles) rather than a FlowSOM code book, so refined leaves
## with heterogeneous marker sets are scored on the same yardstick (§6.1 of
## the design doc: own = distance to the own cluster's representative, rival
## = nearest other representative, ratio > 0.9 = ambiguous).
#'
#' @param X cells x markers (the full-space matrix all arms share)
#' @param R clusters x markers representative matrix
#' @param cl integer cluster id per cell (1..nrow(R))
#' @param idx cell subset to evaluate (memory)
#' @return list(pct_ambiguous, median_ratio, ratio) over idx
amb_cells <- function(X, R, cl, idx = seq_len(nrow(X)), chunk = 50000L) {
    n <- length(idx)
    own <- numeric(n); rival <- numeric(n)
    for (s in seq(1L, n, by = chunk)) {
        i  <- s:min(s + chunk - 1L, n)
        ii <- idx[i]
        Xi <- X[ii, , drop = FALSE]
        D <- outer(rowSums(Xi^2), rowSums(R^2), "+") - 2 * tcrossprod(Xi, R)
        own[i] <- D[cbind(seq_along(i), cl[ii])]
        D[cbind(seq_along(i), cl[ii])] <- Inf
        rival[i] <- rowMins(D)
    }
    ratio <- own / pmax(rival, .Machine$double.eps)
    list(pct_ambiguous = 100 * mean(ratio > 0.9),
         median_ratio = unname(median(ratio)),
         ratio = ratio)
}

#' refine: nearest other cluster (argmin) for each cell of idx. Used by the
#' gap-contribution marker ranking (§5.2 of the design doc).
#' @return integer vector: cluster id of the nearest rival per cell
rival_of <- function(X, R, cl, idx, chunk = 50000L) {
    out <- integer(length(idx))
    for (s in seq(1L, length(idx), by = chunk)) {
        i  <- s:min(s + chunk - 1L, length(idx))
        ii <- idx[i]
        Xi <- X[ii, , drop = FALSE]
        D <- outer(rowSums(Xi^2), rowSums(R^2), "+") - 2 * tcrossprod(Xi, R)
        D[cbind(seq_along(i), cl[ii])] <- Inf
        out[i] <- max.col(-D, ties.method = "first")
    }
    out
}

#' refine: per-marker gap contribution S(m) over a cell subset (§5.2, v1
#' rule; kept for reference — see sep_score for the v2 rule). S(m) =
#' mean_i ( |x_im - c_rival(i),m| - |x_im - c_own,m| ); positive = marker
#' separates the cell from its rival. The marker with the SMALLEST S is the
#' refinement candidate to drop.
#' @param X cells x markers (full space; every marker has a value)
#' @param R clusters x markers representative matrix (cell means)
#' @param cl cluster id per cell
#' @param idx cells to average over (the node's cells)
#' @return named numeric vector, one entry per marker of colnames(X)
gap_scores <- function(X, R, cl, idx) {
    rv <- rival_of(X, R, cl, idx)
    Xs <- X[idx, , drop = FALSE]
    Ro <- R[cl[idx], , drop = FALSE]
    Rr <- R[rv, , drop = FALSE]
    colMeans(abs(Xs - Rr) - abs(Xs - Ro))
}

#' refine: label-free separation score of a k-partition of `Xsub` (v2 WHAT
#' rule, §5.2). Clusters the cells with k-means at k and scores separation as
#' 1 - within_SS / total_SS (between/total variance ratio; ~0 for noise,
#' ~1 for well-separated groups). Recorded 2026-10-28 after the first
#' real-cohort tune run: the v1 smallest-gap rule dropped the marker doing
#' the LEAST rival-separation work, which on real nodes is not the marker
#' carrying the internal mixture (CD19 on T-cell nodes) — the sub-SOM then
#' split a still-mixed cloud into arbitrary halves and 0 of 1692 tuning
#' attempts were accepted. The v2 rule asks directly which marker's REMOVAL
#' best unmakes the mixture, in the space the refinement operates in.
#' @param Xsub cells x markers (the reduced candidate space)
#' @param k number of groups (k_local)
#' @param seed RNG seed for k-means (fixed per candidate for reproducibility)
#' @return numeric in [0, 1]
sep_score <- function(Xsub, k, seed) {
    set.seed(seed)
    km <- suppressWarnings(kmeans(Xsub, centers = k, nstart = 3, iter.max = 25))
    tot <- sum(scale(Xsub, scale = FALSE)^2)
    1 - sum(km$withinss) / tot
}

#' refine: sibling-space ambiguity of a candidate partition (v2 accept
#' criterion, §5.4). The node's cells under `R`/`cl`, measured on the REDUCED
#' marker set only: own/rival among the node's children plus the other leaves
#' restricted to those columns. The refinement must achieve separation in the
#' space the analyst selected; the full-space metric stays REPORTED (its
#' granularity bias is handled by the matched-L control, not by a guard).
#' @param X cells x 12 (full space)
#' @param mk marker columns to score in
#' @param R_full clusters x 12 full-space representative matrix
#' @param cl integer cluster id per cell (full-space partition)
#' @param idx cells to evaluate
sib_amb <- function(X, mk, R_full, cl, idx) {
    amb_cells(X[, mk, drop = FALSE], R_full[, mk, drop = FALSE], cl, idx)
}

## ----------------------------------------------------- preprocessing -------

apply_scaling <- function(X, how) {
    switch(how,
        none          = X,
        center_global = sweep(X, 2, colMeans(X), "-"),
        stop("unknown scaling: ", how))
}

apply_weights <- function(X, how) {
    if (how == "none") return(X)
    if (how == "sqrt_eta2") return(sweep(X, 2, sqrt(eta2(X, pop)), "*"))
    stop("unknown weighting: ", how)
}

prepare_X <- function(Xbase, scaling, weighting) {
    apply_weights(apply_scaling(Xbase, scaling), weighting)
}

#' Nearest SOM node for every cell under distance `distf`, in chunks.
#' Euclidean uses the cross-product identity (only the arg-min is needed).
assign_nodes <- function(X, codes, distf, chunk = 50000L) {
    n <- nrow(X); node <- integer(n)
    for (s in seq(1L, n, by = chunk)) {
        i  <- s:min(s + chunk - 1L, n)
        Xi <- X[i, , drop = FALSE]
        D <- if (distf == 2) {
            outer(rowSums(Xi^2), rowSums(codes^2), "+") - 2 * tcrossprod(Xi, codes)
        } else node_distances(Xi, codes, distf)
        node[i] <- max.col(-D, ties.method = "first")
    }
    node
}

## ---------------------------------------------------------- run engine ----

cfg_cols  <- c("cofactor", "batch", "scaling", "weighting", "distf", "markers")
run_cols  <- c(cfg_cols, "seed", "xdim", "k")
n_seeds   <- 5L
seeds_all <- 1234L + seq_len(n_seeds) - 1L

run_key <- function(g) {
    gsub("[^A-Za-z0-9_.-]", "_", paste(g$cofactor, g$batch, g$scaling, g$weighting,
        paste0("d", g$distf), paste0("m", g$markers), paste0("s", g$seed),
        paste0("g", g$xdim), paste0("k", g$k), sep = "__"))
}

#' Run one configuration (flat SOM; the detailed study's run_one with the
#' marker set made explicit and the seed-controlled mapping-agreement check
#' kept).
#' @param Xbase cells x type-markers after cofactor (already subset to
#'   cfg$markers)
#' @param cfg one-row data.frame with the columns in `run_cols`
#' @return one-row data.frame: configuration + metrics.
run_one <- function(Xbase, cfg) {
    t0 <- proc.time()[["elapsed"]]
    out <- data.frame(cfg[, run_cols], ari = NA_real_, ari_node_majority = NA_real_,
        ari_donor = NA_real_, node_purity = NA_real_, pct_ambiguous = NA_real_,
        median_ratio = NA_real_, pct_mapping_agree = NA_real_, n_empty_nodes = NA_integer_,
        runtime_s = NA_real_, error = NA_character_, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    tryCatch({
        X <- prepare_X(Xbase, cfg$scaling, cfg$weighting)
        set.seed(cfg$seed)
        som <- FlowSOM::SOM(X, xdim = cfg$xdim, ydim = cfg$xdim, rlen = 10, mst = 1,
                            distf = cfg$distf, silent = TRUE)
        stopifnot(!is.null(som$mapping))
        node <- assign_nodes(X, som$codes, cfg$distf)
        meta <- meta_cluster(som$codes, cfg$k, cfg$seed)
        cl   <- meta[node]
        maj  <- apply(table(node, pop), 1, which.max)
        out$ari               <- ari_score(cl, pop)
        out$ari_node_majority <- ari_score(maj[as.character(node)], as.integer(pop))
        out$ari_donor         <- ari_score(cl, donor)
        out$node_purity       <- node_purity(node, pop)
        out$n_empty_nodes     <- nrow(som$codes) - length(unique(node))
        amb <- ambiguity_metrics(X, som$codes, node, cfg$distf, amb_idx)
        out$pct_ambiguous     <- amb$pct_ambiguous
        out$median_ratio      <- amb$median_ratio
        out$pct_mapping_agree <- 100 * mean(as.integer(som$mapping[, 1]) == node)
    }, error = function(e) out$error <<- conditionMessage(e))
    out$runtime_s <- proc.time()[["elapsed"]] - t0
    out
}

#' refine: one flat-SOM run that ALSO returns the partition and cell-mean
#' representatives, for the mean-based yardstick (§6.1).
#' @return list(metrics = one-row data.frame like run_one's, cl = integer
#'   cluster id per cell (metaclusters), reps = k x ncol(X) cell-mean matrix)
flat_mean_run <- function(Xbase, cfg) {
    X <- prepare_X(Xbase, cfg$scaling, cfg$weighting)
    t0 <- proc.time()[["elapsed"]]
    set.seed(cfg$seed)
    som <- FlowSOM::SOM(X, xdim = cfg$xdim, ydim = cfg$xdim, rlen = 10, mst = 1,
                        distf = cfg$distf, silent = TRUE)
    node <- assign_nodes(X, som$codes, cfg$distf)
    meta <- meta_cluster(som$codes, cfg$k, cfg$seed)
    cl   <- meta[node]
    maj  <- apply(table(node, pop), 1, which.max)
    amb  <- ambiguity_metrics(X, som$codes, node, cfg$distf, amb_idx)
    ## Mean-based yardstick on the same partition (§6.1): reps are per-cluster
    ## rows (amb_cells maps cells to them via cl). Meta-cluster ids may skip
    ## values when a cluster consists only of unpopulated SOM codes (observed
    ## at k = 102: 101 non-empty clusters among ids 1..102) — renumber the
    ## observed ids to contiguous 1..n_distinct so the reps rows and cl agree.
    ## The partition itself is unchanged, only its labels; ARI, purity and
    ## the mean-based ambiguity are label-invariant.
    cl_int <- as.integer(cl)
    if (anyNA(cl_int) || min(cl_int) < 1L ||
        length(unique(cl_int)) < max(cl_int)) {
        cl_int <- match(cl_int, sort(unique(cl_int)))
    }
    reps <- do.call(rbind, lapply(split(seq_len(N), cl_int),
                                  function(ii) colMeans(X[ii, , drop = FALSE])))
    m <- data.frame(cfg[, run_cols], ari = ari_score(cl, pop),
        ari_node_majority = ari_score(maj[as.character(node)], as.integer(pop)),
        ari_donor = ari_score(cl, donor), node_purity = node_purity(node, pop),
        pct_ambiguous = amb$pct_ambiguous, median_ratio = amb$median_ratio,
        pct_mapping_agree = 100 * mean(as.integer(som$mapping[, 1]) == node),
        n_empty_nodes = nrow(som$codes) - length(unique(node)),
        runtime_s = proc.time()[["elapsed"]] - t0,
        stringsAsFactors = FALSE)
    ## mean-based yardstick on the same partition (§6.1): reps and cl_int
    ## agree after the renumbering above
    amb_mean <- amb_cells(X, reps, cl_int, amb_idx)
    m$pct_ambiguous_mean <- amb_mean$pct_ambiguous
    m$median_ratio_mean  <- amb_mean$median_ratio
    list(metrics = m, cl = cl_int, reps = reps)
}

#' Run (or resume) a stage of flat runs.
run_stage <- function(grid, stage) {
    stopifnot(all(run_cols %in% names(grid)))
    dir <- file.path(out_dir, "runs", stage)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    grid$file <- file.path(dir, paste0(run_key(grid), ".rds"))
    todo <- which(!file.exists(grid$file))
    message(sprintf("[%s] %d runs: %d cached, %d to compute (%d core(s))",
        stage, nrow(grid), nrow(grid) - length(todo), length(todo), n_cores))
    groups <- split(todo, paste(grid$cofactor[todo], grid$batch[todo],
                                grid$markers[todo], sep = "||"))
    for (g in names(groups)) {
        rows <- groups[[g]]
        message(sprintf("[%s] base matrix '%s': %d runs", stage, g, length(rows)))
        Xbase <- get_base(grid$cofactor[rows[1]], grid$markers[rows[1]])
        worker <- function(i) {
            res <- run_one(Xbase, grid[i, ])
            save_rds_atomic(res, grid$file[i])
            invisible(NULL)
        }
        invisible(parallel::mclapply(rows, worker, mc.cores = n_cores, mc.preschedule = FALSE))
        rm(Xbase); invisible(gc())
    }
    res <- lapply(grid$file, function(f) if (file.exists(f)) readRDS(f) else NULL)
    if (any(vapply(res, is.null, NA))) {
        warning(sum(vapply(res, is.null, NA)), " run(s) produced no result file ",
                "(worker died, e.g. out of memory). Re-run the chunk to retry them.")
    }
    res <- do.call(rbind, res)
    res$distance <- factor(dist_labels[as.character(res$distf)],
                           levels = unname(dist_labels)[c(2, 1, 3, 4)])
    rownames(res) <- NULL
    res
}

## ------------------------------------------------------ refinement tree ----

#' refine: default refinement parameters (§5, §7 of the design doc). The
#' freeze protocol tunes tau_refine/tau_stop/min_cells on seed 1234 only.
#' v3 semantics (2026-10-28, recorded in §11.2 of the design doc):
#' accept_dpur = MINIMUM purity gain for acceptance (was: maximum purity
#' degradation), accept_eps = MAXIMUM tolerated worsening of the node's
#' full-space ambiguity (was: minimum required improvement). Rationale: the
#' mean-based own/rival ambiguity metric is structurally anti-refinement —
#' splitting a node adds the sibling as the closest rival, and populations
#' inside one SOM node are by definition the closest to each other, so
#' ambiguity can almost never improve under splitting (0 of 1692 v1 tune
#' attempts and 0 of 9 v2 sibling-space attempts were accepted; purity gains
#' of +0.09..+0.15 from genuine unmixing were invisible to both rules). The
#' refinement's goal is therefore measured in purity (the WHERE premise's
#' own diagnostic), with the geometry cost bounded.
refine_params <- function(tau_refine = 15, purity_refine = 0.9, tau_stop = 10,
                          min_cells = 2000, min_markers = 6, max_depth = 3,
                          accept_dpur = 0.02, accept_eps = 1) {
    list(tau_refine = tau_refine, purity_refine = purity_refine,
         tau_stop = tau_stop, min_cells = min_cells, min_markers = min_markers,
         max_depth = max_depth, accept_dpur = accept_dpur, accept_eps = accept_eps)
}

#' refine: build the hierarchical refinement tree for one seed.
#'
#' Level 0 = the anchor flat SOM (10x10, k=12). Nodes qualify for refinement
#' by CySA's own diagnostics (leaf ambiguity > tau_refine OR purity <
#' purity_refine) and size (>= min_cells). Each attempt drops one marker
#' (v2: the one whose removal maximizes k-means separation of the node's
#' cells in the reduced space; "gap" = the v1 rule kept for reference,
#' "random" = the ablation), trains a sub-SOM on the node's cells in the
#' reduced space, metaclusters it into k_local groups, and accepts iff the
#' node's purity improves by at least accept_dpur (v3: unmixing achieved)
#' and its full-space ambiguity worsens by at most accept_eps (bounded
#' geometry cost). Sibling-space ambiguity is logged before/after as the
#' reported cost of splitting (§11.2), and the cohort-wide ambiguity too;
#' neither is gated. Rejected attempts keep the parent as a final leaf.
#'
#' @param Xbase cells x 12 prepared matrix (the full-space yardstick space)
#' @param seed RNG seed for the level-0 SOM
#' @param p refinement parameters (refine_params())
#' @param marker_choice "sep" (v2, default), "gap" (v1, kept for reference),
#'   or "random" (ablation arm)
#' @return list(leaves, leaf_of, log, metrics, R, runtime_s)
refine_tree <- function(Xbase, seed, p = refine_params(), marker_choice = "sep") {
    t0 <- proc.time()[["elapsed"]]
    stopifnot(marker_choice %in% c("sep", "gap", "random"))
    X <- Xbase
    markers_all <- colnames(X)
    N_loc <- nrow(X)

    ## level 0: anchor flat SOM
    set.seed(seed)
    som0 <- FlowSOM::SOM(X, xdim = 10, ydim = 10, rlen = 10, mst = 1, distf = 2,
                         silent = TRUE)
    node0 <- assign_nodes(X, som0$codes, 2)
    meta0 <- meta_cluster(som0$codes, 12, seed)
    occ <- sort(unique(node0))
    leaves <- lapply(occ, function(nd) list(
        idx = which(node0 == nd), markers = markers_all, depth = 0L,
        parent = NA_integer_, terminal = FALSE, child_of = nd))
    leaf_of <- integer(N_loc); for (l in seq_along(leaves)) leaf_of[leaves[[l]]$idx] <- l
    for (l in seq_along(leaves)) leaves[[l]]$mean <- colMeans(X[leaves[[l]]$idx, , drop = FALSE])

    rep_mat <- function(ls) do.call(rbind, lapply(ls, `[[`, "mean"))
    log <- list()
    child_seed_base <- seed * 1000L
    n_children_made <- 0L

    ## NOTE: regular <- never modifies a binding in an enclosing frame, so a
    ## nested function must take and RETURN the leaf list explicitly instead
    ## of mutating it in place.
    refresh_stats <- function(ls) {
        R <- rep_mat(ls)
        amb <- amb_cells(X, R, leaf_of)
        for (l in seq_along(ls)) {
            ii <- ls[[l]]$idx
            lv <- ls[[l]]
            lv$pct_amb  <- 100 * mean(amb$ratio[ii] > 0.9)
            lv$purity   <- node_purity_one(ii)
            lv$median_r <- unname(median(amb$ratio[ii]))
            ls[[l]] <- lv
        }
        ls
    }
    leaves <- refresh_stats(leaves)

    ## Stopping: tau_stop floors the AMBIGUITY branch only — a node that is
    ## already confident (pct_amb < tau_stop) is terminal unless it is impure
    ## (purity < purity_refine keeps refinement alive for confidently-mixed
    ## nodes). With tau_stop <= tau_refine the condition is a no-op, so those
    ## parameter combos behave exactly as before.
    refinable <- function(l) {
        lv <- leaves[[l]]
        !isTRUE(lv$terminal) && length(lv$idx) >= p$min_cells &&
            lv$depth < p$max_depth && length(lv$markers) > p$min_markers &&
            (lv$pct_amb > p$tau_refine || lv$purity < p$purity_refine) &&
            (lv$pct_amb >= p$tau_stop || lv$purity < p$purity_refine)
    }
    queue <- which(vapply(seq_along(leaves), refinable, NA))
    queue_key <- vapply(leaves[queue], function(lv) lv$pct_amb, 0)

    grid_side <- function(n) {
        nodes <- max(9, round(n / (N_loc / 100)))
        min(10, max(3, ceiling(sqrt(nodes))))
    }

    while (length(queue)) {
        top <- which.max(queue_key)
        l  <- queue[top]; queue <- queue[-top]; queue_key <- queue_key[-top]
        lv <- leaves[[l]]
        if (!refinable(l)) next
        ii  <- lv$idx
        n_i <- length(ii)
        R_now <- rep_mat(leaves)
        cand  <- lv$markers
        ## full-space before (reported; the cross-arm yardstick)
        before <- amb_cells(X, R_now, leaf_of, ii)
        pur_before <- node_purity_one(ii)
        k_tab <- table(pop_int[ii])
        k_local <- min(sum(k_tab >= 0.05 * sum(k_tab)), 8L)
        child_seed <- child_seed_base + 1000L * (lv$depth + 1L) + length(log) + 1L
        gs <- grid_side(n_i)
        row <- data.frame(step = length(log) + 1L, seed = seed, choice = marker_choice,
            depth = lv$depth, parent = l, n_cells = n_i, marker_dropped = NA_character_,
            score_drop = NA_real_, k_local = k_local, grid = gs,
            pct_amb_before = before$pct_ambiguous, pur_before = pur_before,
            pct_amb_after = NA_real_, pur_after = NA_real_,
            pct_amb_sib_before = NA_real_, pct_amb_sib_after = NA_real_,
            pct_amb_after_local = NA_real_,
            pct_amb_global_before = NA_real_, pct_amb_global_after = NA_real_,
            accepted = NA, stringsAsFactors = FALSE)
        if (k_local < 2L) {
            row$accepted <- FALSE
            leaves[[l]]$terminal <- TRUE
            log[[length(log) + 1L]] <- row
            next
        }

        ## choose the marker to drop (§5.2; v2 rule: separation-maximizing)
        m_drop <- NA_character_; score_drop <- NA_real_
        if (marker_choice == "sep") {
            scores <- vapply(seq_along(cand), function(q)
                sep_score(X[ii, setdiff(cand, cand[q]), drop = FALSE], k_local,
                          child_seed + 5000L + q), 0)
            names(scores) <- cand
            m_drop <- names(which.max(scores))
            score_drop <- unname(scores[m_drop])
        } else if (marker_choice == "gap") {
            S <- gap_scores(X, R_now, leaf_of, ii)[cand]
            m_drop <- names(which.min(S))
            score_drop <- unname(S[m_drop])
        } else {
            set.seed(child_seed + 7000L + length(log))
            m_drop <- sample(cand, 1)
        }

        mk <- setdiff(cand, m_drop)
        Xsub <- X[ii, mk, drop = FALSE]
        row$marker_dropped <- m_drop
        row$score_drop <- score_drop
        set.seed(child_seed)
        som <- FlowSOM::SOM(Xsub, xdim = gs, ydim = gs, rlen = 10, mst = 1, distf = 2,
                            silent = TRUE)
        node_sub <- assign_nodes(Xsub, som$codes, 2)
        meta <- meta_cluster(som$codes, k_local, child_seed)
        cl_sub <- meta[node_sub]
        ## code-based sibling-space ambiguity of the sub-SOM (§6.3)
        local_amb <- ambiguity_metrics(Xsub, som$codes, node_sub, 2, seq_len(n_i))

        parts <- split(seq_len(n_i), cl_sub)
        ## candidate partition: leaf l is replaced by its children. New-id
        ## space: leaves 1..l-1 keep ids; leaves l+1..L shift down by one;
        ## children take L-1+1 .. L-1+q. R_after rows follow the same order.
        L_now <- length(leaves)
        R_after <- rbind(R_now[-l, , drop = FALSE],
                         do.call(rbind, lapply(parts, function(jj) colMeans(X[ii[jj], , drop = FALSE]))))
        cl_after <- leaf_of
        if (l < L_now) cl_after[cl_after > l] <- cl_after[cl_after > l] - 1L
        for (q in seq_along(parts)) cl_after[ii[parts[[q]]]] <- L_now - 1L + q
        after <- amb_cells(X, R_after, cl_after, ii)
        pur_after <- mean(vapply(parts, function(jj) node_purity_one(ii[jj]), 0))
        ## full-space cohort-wide effect: REPORTED, not gated (the matched-L
        ## control A3 handles the granularity bias; §5.4 v2)
        amb_g_before <- amb_cells(X, R_now, leaf_of, amb_idx)
        amb_g_after  <- amb_cells(X, R_after, cl_after, amb_idx)
        ## sibling-space before/after (REPORTED diagnostic: the geometry cost
        ## of splitting — the sibling becomes the closest rival, §11.2)
        before_sib <- amb_cells(X[, cand, drop = FALSE], R_now[, cand, drop = FALSE],
                                leaf_of, ii)
        after_sib  <- amb_cells(X[, mk, drop = FALSE], R_after[, mk, drop = FALSE],
                                cl_after, ii)

        row$pct_amb_after <- after$pct_ambiguous
        row$pur_after <- pur_after
        row$pct_amb_sib_before <- before_sib$pct_ambiguous
        row$pct_amb_sib_after  <- after_sib$pct_ambiguous
        row$pct_amb_after_local <- local_amb$pct_ambiguous
        row$pct_amb_global_before <- amb_g_before$pct_ambiguous
        row$pct_amb_global_after  <- amb_g_after$pct_ambiguous
        ## v3 accept (§5.4 v3): purity must GAIN by at least accept_dpur
        ## (unmixing achieved); the node's full-space ambiguity may worsen by
        ## at most accept_eps (bounded geometry cost). Sibling-space columns
        ## above are reported diagnostics (the sibling-rivalry cost), not gates.
        ok <- pur_after >= pur_before + p$accept_dpur &&
              after$pct_ambiguous <= before$pct_ambiguous + p$accept_eps
        row$accepted <- ok
        log[[length(log) + 1L]] <- row
        if (!ok) {
            leaves[[l]]$terminal <- TRUE
            next
        }
        n_children_made <- n_children_made + length(parts)
        new_leaves <- lapply(seq_along(parts), function(q) list(
            idx = ii[parts[[q]]], markers = mk, depth = lv$depth + 1L,
            parent = l, terminal = FALSE, child_of = q))
        ## commit: same id mapping as the evaluation above — children take the
        ## END ids (L_now-1+q, matching cl_after), so the leaves list must
        ## append them at the end, not insert at position l. Means are set on
        ## new_leaves BEFORE the list is rebuilt so `leaves` shares them.
        for (q in seq_along(new_leaves)) new_leaves[[q]]$mean <- R_after[L_now - 1L + q, ]
        leaves <- c(leaves[-l], new_leaves)
        leaf_of <- cl_after
        leaves <- refresh_stats(leaves)
        queue[queue > l] <- queue[queue > l] - 1L
        for (q in seq_along(new_leaves)) {
            lid <- L_now - 1L + q
            if (refinable(lid)) { queue <- c(queue, lid); queue_key <- c(queue_key, leaves[[lid]]$pct_amb) }
        }
    }

    R <- rep_mat(leaves)
    cl <- leaf_of
    amb <- amb_cells(X, R, cl, amb_idx)
    ## majority LABEL per leaf; ari_node_majority is the anchor-comparable
    ## ARI family member (the anchor's 0.9872 is computed the same way)
    maj_lab <- vapply(seq_along(leaves), function(l2)
        names(which.max(table(pop[leaves[[l2]]$idx]))), "")
    pur_total <- sum(vapply(leaves, function(l2)
        max(tabulate(pop_int[l2$idx])), 0L)) / N_loc
    metrics <- data.frame(seed = seed, choice = marker_choice,
        tau_refine = p$tau_refine, tau_stop = p$tau_stop, min_cells = p$min_cells,
        accept_eps = p$accept_eps,
        L = length(leaves), ari = ari_score(cl, pop),
        ari_node_majority = ari_score(maj_lab[cl], as.character(pop)),
        node_purity = pur_total,
        pct_ambiguous = amb$pct_ambiguous, median_ratio = amb$median_ratio,
        n_attempts = length(log),
        n_accepted = sum(vapply(log, function(r) isTRUE(r$accepted), NA)),
        runtime_s = proc.time()[["elapsed"]] - t0,
        stringsAsFactors = FALSE)
    list(leaves = leaves, leaf_of = leaf_of, log = do.call(rbind, log),
         metrics = metrics, R = R, runtime_s = metrics$runtime_s)
}

#' refine: expert-gating upper bound (A4). Representatives = population cell
#' means in the full space; purity/ARI trivially 1 (the partition IS the
#' labels). The ambiguity metric is the interesting part.
expert_bound <- function(X) {
    pops <- sort(unique(pop))
    reps <- do.call(rbind, lapply(pops, function(pp) colMeans(X[pop == pp, , drop = FALSE])))
    cl <- as.integer(pop)[match(pop, pops)]
    amb <- amb_cells(X, reps, cl, amb_idx)
    data.frame(partition = "expert_gating", L = length(pops), ari = 1,
               node_purity = 1, pct_ambiguous = amb$pct_ambiguous,
               median_ratio = amb$median_ratio, stringsAsFactors = FALSE)
}

## --------------------------------------------------------- stage drivers ---

## Stage drivers run on the FULL cohort only (dev renders read caches; a
## fresh dev-mode run would poison the full-data cache keys).
dev_guard <- function(what) {
    if (dev_n > 0) stop("stage driver '", what,
                        "' runs on the full cohort only; dev renders read caches")
}

#' S1/A1: 12-marker anchor runs (the d5-7 CCR7 exclusion applied to the
#' anchor config), code-based AND mean-based metrics per seed.
stage_anchor12 <- function(seeds = seeds_all, xdim = 10, k = 12) {
    dev_guard("anchor12")
    Xbase <- cached("00_base_paper12.rds", get_base("paper", type_m12), csv = FALSE)
    dir.create(file.path(out_dir, "anchor"), recursive = TRUE, showWarnings = FALSE)
    out <- lapply(seeds, function(sv) {
        cfg <- data.frame(cofactor = "paper", batch = "none", scaling = "none",
                          weighting = "none", distf = 2, markers = "type_m12",
                          seed = sv, xdim = xdim, k = k, stringsAsFactors = FALSE)
        f <- file.path(out_dir, "anchor", paste0(run_key(cfg), ".rds"))
        if (!file.exists(f)) save_rds_atomic(flat_mean_run(Xbase, cfg), f)
        readRDS(f)
    })
    do.call(rbind, lapply(out, `[[`, "metrics"))
}

#' S1: the 13-marker side of the paired comparison (CCR7 included, i.e. the
#' d5-7 artifact NOT excluded), same anchor config and protocol as
#' stage_anchor12 so the pairing is like-for-like. Cached under
#' markers = "type_m" keys, separate from the 12-marker anchors.
stage_anchor13 <- function(seeds = seeds_all, xdim = 10, k = 12) {
    dev_guard("anchor13")
    Xbase <- cached("00_base_paper13.rds", get_base("paper", type_m), csv = FALSE)
    dir.create(file.path(out_dir, "anchor"), recursive = TRUE, showWarnings = FALSE)
    out <- lapply(seeds, function(sv) {
        cfg <- data.frame(cofactor = "paper", batch = "none", scaling = "none",
                          weighting = "none", distf = 2, markers = "type_m",
                          seed = sv, xdim = xdim, k = k, stringsAsFactors = FALSE)
        f <- file.path(out_dir, "anchor", paste0(run_key(cfg), ".rds"))
        if (!file.exists(f)) save_rds_atomic(flat_mean_run(Xbase, cfg), f)
        readRDS(f)
    })
    do.call(rbind, lapply(out, `[[`, "metrics"))
}

#' Identity baseline: the level-0 anchor's populated SOM codes with per-code
#' mean representatives — what the refinement tree starts from before any
#' refinement. Scored with the same amb_cells family as the tree metrics
#' (§6.1), per seed, so the arms table can compare A2 against its own
#' un-refined start. The SOM is rebuilt per seed exactly as refine_tree()'s
#' level-0 does (same seed -> same codes -> same partition).
stage_identity <- function(seeds = seeds_all) {
    dev_guard("identity")
    Xbase <- cached("00_base_paper12.rds", get_base("paper", type_m12), csv = FALSE)
    X <- prepare_X(Xbase, "none", "none")
    out <- do.call(rbind, lapply(seeds, function(sv) {
        set.seed(sv)
        som0 <- FlowSOM::SOM(X, xdim = 10, ydim = 10, rlen = 10, mst = 1,
                             distf = 2, silent = TRUE)
        node0 <- assign_nodes(X, som0$codes, 2)
        occ <- sort(unique(node0))
        leaf_of <- integer(nrow(X))
        for (l in seq_along(occ)) leaf_of[node0 == occ[l]] <- l
        R <- do.call(rbind, lapply(occ, function(nd)
            colMeans(X[node0 == nd, , drop = FALSE])))
        amb <- amb_cells(X, R, leaf_of)
        pur <- sum(vapply(occ, function(nd)
            max(tabulate(pop_int[node0 == nd])), 0L)) / nrow(X)
        maj <- vapply(seq_along(occ), function(l)
            names(which.max(table(pop[node0 == occ[l]]))), "")
        ## column names follow the arms-table vocabulary (§6.4) so the
        ## cached file and this regeneration path stay interchangeable
        data.frame(seed = sv, L = length(occ), ari = ari_score(leaf_of, pop),
                   ari_node_majority = ari_score(maj[leaf_of], as.character(pop)),
                   purity = pur, pct_amb = amb$pct_ambiguous,
                   median_ratio = amb$median_ratio, stringsAsFactors = FALSE)
    }))
    save_rds_atomic(out, file.path(out_dir, "identity_table.rds"))
    message("[STATS] identity rows: ", nrow(out))
    out
}

#' Parameter freeze protocol (§7): trees on seed 1234 over a small threshold
#' grid. Returns one row per combination; the frozen values are chosen from
#' it and recorded in the vignette.
stage_tune <- function(combos, seed = 1234L) {
    dev_guard("tune")
    Xbase <- cached("00_base_paper12.rds", get_base("paper", type_m12), csv = FALSE)
    X <- prepare_X(Xbase, "none", "none")
    combo_key <- function(r) paste(sprintf("tune_v3_s%d", seed),
        paste(paste(names(combos), as.character(unlist(combos[r, ])),
                    sep = "_"), collapse = "_"),
        sep = "_")
    one <- function(r) {
        p <- do.call(refine_params, as.list(combos[r, ]))
        key <- combo_key(r)
        f <- file.path(out_dir, "runs", paste0(key, ".rds"))
        if (!file.exists(f)) {
            res <- refine_tree(X, seed, p, "sep")
            save_rds_atomic(res, f)
        } else res <- readRDS(f)
        res$metrics
    }
    res_l <- parallel::mclapply(seq_len(nrow(combos)), one, mc.cores = n_cores,
                                mc.preschedule = FALSE)
    crashed <- which(vapply(res_l, is.null, TRUE))
    if (length(crashed)) stop("stage_tune: ", length(crashed),
                              " combo(s) crashed in workers; see the job log")
    out <- do.call(rbind, res_l)
    rownames(out) <- NULL
    out
}

#' A2 + A5: refinement trees, guided and random, all seeds.
stage_trees <- function(seeds = seeds_all, p = refine_params()) {
    dev_guard("trees")
    Xbase <- cached("00_base_paper12.rds", get_base("paper", type_m12), csv = FALSE)
    X <- prepare_X(Xbase, "none", "none")
    dir.create(file.path(out_dir, "trees"), recursive = TRUE, showWarnings = FALSE)
    worker <- function(sv) {
        for (choice in c("sep", "random")) {
            f <- file.path(out_dir, "trees",
                           sprintf("tree_v3_s%d_%s.rds", sv, choice))
            if (file.exists(f)) next
            res <- refine_tree(X, sv, p, choice)
            save_rds_atomic(res, f)
        }
        invisible(NULL)
    }
    invisible(parallel::mclapply(seeds, worker, mc.cores = n_cores, mc.preschedule = FALSE))
    do.call(rbind, lapply(seeds, function(sv) {
        rbind(readRDS(file.path(out_dir, "trees", sprintf("tree_v3_s%d_sep.rds", sv)))$metrics,
              readRDS(file.path(out_dir, "trees", sprintf("tree_v3_s%d_random.rds", sv)))$metrics)
    }))
}

#' A3: flat matched controls. For each seed: flat 20x20 metaclustered at
#' k = L_tree (partition-level match) and flat 20x20 at k = 12 (same number
#' of final clusters as the anchor, more nodes).
stage_controls <- function(seeds = seeds_all) {
    dev_guard("controls")
    Xbase <- cached("00_base_paper12.rds", get_base("paper", type_m12), csv = FALSE)
    for (sv in seeds) {
        tree <- readRDS(file.path(out_dir, "trees", sprintf("tree_v3_s%d_sep.rds", sv)))
        L_tree <- tree$metrics$L
        grid <- merge(
            data.frame(cofactor = "paper", batch = "none", scaling = "none",
                       weighting = "none", distf = 2, markers = "type_m12",
                       xdim = 20, k = c(L_tree, 12), stringsAsFactors = FALSE),
            data.frame(seed = sv), by = NULL)
        for (r in seq_len(nrow(grid))) {
            key <- run_key(grid[r, ])
            f <- file.path(out_dir, "runs", "controls", paste0(key, ".rds"))
            dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
            if (!file.exists(f)) save_rds_atomic(flat_mean_run(Xbase, grid[r, ]), f)
        }
    }
    invisible(NULL)
}

#' A4: expert-gating upper bound, cached (the vignette reads it; the SLURM
#' driver computes it once on the full cohort).
stage_expert <- function() {
    dev_guard("expert")
    Xbase <- cached("00_base_paper12.rds", get_base("paper", type_m12), csv = FALSE)
    X <- prepare_X(Xbase, "none", "none")
    cached("a4_expert_bound.rds", expert_bound(X))
}
