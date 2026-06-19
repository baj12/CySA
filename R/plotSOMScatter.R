#' @rdname plotScatter
#' @title Scatter plot
#' 
#' @description Bivariate scatter plots including visualization of
#' (group-specific) gates, their boundaries and percentage of selected cells.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character vector specifying which channels to plot.
#' @param metaSlot name of the metadata slot containing SOM codes.
#' @param pointSize column in SOM stats used to size points.
#' @param color_by column to color points by.
#' @param bins number of bins when coloring by density.
#' @param assay name of the assay to use.
#' @param statsSlot name of the metadata slot containing SOM stats.
#' @param label how axis labels should be constructed.
#' @param zeros logical specifying whether to include 0 values.
#' @param k_pal optional cluster color palette.
#' @param xRN optional row names.
#' @param xCN optional channel names.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' \dontrun{
#' # Requires a SingleCellExperiment with SOM_codes and SOM_stats metadata.
#' plotSOMScatter(sce, chs = c("FSC.A", "SSC.A"))
#' }
#'
#' @export
plotSOMScatter <- function(x, chs, metaSlot = "SOM_codes", pointSize = "n",
                           color_by = "n",
                           bins = 100, assay = "exprs", statsSlot="SOM_stats",
                           label = c("target", "channel", "both"),
                           zeros = FALSE, k_pal = NULL,
                           xRN = NULL, xCN=NULL) {
  # check validity of input arguments
  label <- match.arg(label)
  # Ensure pointSize and color_by are single values
  pointSize <- pointSize[1]
  color_by <- color_by[1]
  args <- as.list(environment())
  # CATALYST:::.check_args_plotScatter(args)
  if(!metaSlot %in% names(S4Vectors::metadata(x))) {
    stop("Need ", metaSlot, " in metadata of sce")
  }
  if(!statsSlot %in% names(S4Vectors::metadata(x))) {
    warning(statsSlot, " not found in metadata of sce - proceeding without stats")
    statsSlot <- NULL
  }
  # 2do apply parallel
  # compute stats if requested, this is time consuming
  for(ch in chs){
    if(!ch %in% colnames(S4Vectors::metadata(x)[[metaSlot]])){
      warning("computing stats for ", ch, "\n")
      mdt =S4Vectors::metadata(x)[[metaSlot]]
      # mdt[,ch] = 0
      mdt = cbind(mdt, matrix(0, nrow = nrow(mdt), ncol = 1))
      colnames(mdt)[ncol(mdt)] <- ch
      mdtt =S4Vectors::metadata(x)
      for(cl in as.integer(unique(colData(x)$cluster_id))){
        mdt[cl,ch] = mean(assays(x)[[assay]][ch,which(colData(x)$cluster_id == cl)])
      }
      mdtt[[metaSlot]]= mdt
      S4Vectors::metadata(x) = mdtt
    }
  }
  # subset features to speed up matrix transpose  ‘
  if(is.null(xRN))
    xRN <- rownames(x)
  if(is.null(xCN))
    xCN <- CATALYST::channels(x)
  i <- lapply(list(xRN, xCN), function(u) {
    i <- match(chs, u, nomatch = 0)
    if (all(i == 0)) NULL else i
  })
  i <- unlist(i)

  # Check if requested channels exist in SOM_codes
  som_codes <- S4Vectors::metadata(x)[[metaSlot]]
  missing_chs <- setdiff(chs, colnames(som_codes))
  if (length(missing_chs) > 0) {
    warning("Channels not found in SOM_codes: ", paste(missing_chs, collapse = ", "),
            "\nAvailable: ", paste(colnames(som_codes), collapse = ", "))
    # Filter to only available channels
    chs <- intersect(chs, colnames(som_codes))
    if (length(chs) < 2) {
      return(NULL)  # Need at least 2 channels for scatter plot
    }
  }
  yy <- som_codes[, chs, drop = FALSE]
  # y <- x[i, , drop = FALSE]
  # y <- x[chs, , drop = FALSE] # this seems to be an expensive operation
  # y <- assay(x, assay)
  # y = y[chs, , drop = FALSE]
  # rename features for visualization to
  # include both channel name & description
  nms <- switch(label, target = xRN, channel = xCN,
                both = ifelse(xCN == xRN, xCN, paste(xCN, xRN, sep = "-")))
  # chs[i != 0] <- rownames(y) <- nms[i]
  
  # construct data.frame of specified assay data & all cell metadata
  if (isTRUE(color_by %in% names(CATALYST::cluster_codes(x))))
    x[[color_by]] <- CATALYST::cluster_ids(x, color_by)
  cd <- cbind(SummarizedExperiment::colData(x), SingleCellExperiment::int_colData(x))


  # df <- data.frame(
  #   t(as.matrix(y)), cd,
  #   check.names = FALSE,
  #   stringsAsFactors = FALSE)
  # cd_vars <- intersect(names(cd), names(df))


  if (!is.null(statsSlot)) {
    stats <- S4Vectors::metadata(x)[[statsSlot]]
    df <- cbind(yy, stats)
  } else {
    df <- yy
    # Add a dummy n column for pointSize when stats are missing
    df$n <- 1
  }
  # qualify rowSums to avoid ambiguity
  if (!zeros) df <- df[Matrix::rowSums(df[, chs[c(1, 2)]] == 0) == 0, ]
  # browser()
  # initialize faceting & (optionally) melt data.frame
  if (length(chs) > 2) {
    # df <- melt(df, id.vars = unique(c(chs[1], cd_vars)))
    if (!is.null(statsSlot)) {
      df <- reshape2::melt(df, id.vars = unique(c(chs[1], colnames(stats))))
    } else {
      df <- reshape2::melt(df, id.vars = chs[1])
    }
    facet <- "variable"
    ylab <- ylab(NULL)
    chs[2] <- "value"
  } else {
    facet <- NULL
    ylab <- NULL
  }
  
  fill_var <- NULL
  col_var <- sprintf("%s", color_by)

  # If statsSlot is NULL, check if color_by column exists in df
  if (!is.null(statsSlot) || !col_var %in% colnames(df)) {
    # Default to density coloring if stats not available
    col_var <- NULL
    fill_var <- NULL
  }

  geom <- geom_point(alpha = 0.2, na.rm = TRUE)
  guides <- NULL
  scales <- NULL

  if (!is.null(col_var) && col_var %in% colnames(df) && is.numeric(df[[col_var]])) {
    # Continuous color scale - don't set guides to avoid conflict
    scales <- scale_color_gradientn(
      colors = c("navy", rev(RColorBrewer::brewer.pal(11, "Spectral"))),
      guide = "colorbar"
    )
  } else if (!is.null(col_var) && col_var %in% names(CATALYST::cluster_codes(x))) {
    if (is.null(k_pal)) k_pal <- CySA_default_cluster_cols()
    scales <- scale_color_manual(values = k_pal)
    guides <- guides(col = guide_legend(
      override.aes = list(alpha = 1, size = 3)))
  }
  
  # facet <- c(facet, facet_by)
  # if (!is.null(facet)) {
  #   if (length(facet) == 1) {
  #     facet <- facet_wrap(facet)    
  #   } else {
  #     facet <- facet_grid(
  #       cols = vars(!!sym(facet[1])), 
  #       rows = vars(!!sym(facet[2])))
  #   }
  # }
  if(is.null(chs))return(NULL)
  xy <- sprintf("%s", chs)
  if (!zeros) df <- df[Matrix::rowSums(df[, chs[c(1, 2)]] == 0) == 0, ]

  # Build aesthetic mapping using aes_string for proper evaluation
  aes_list <- aes_string(
    x = xy[1],
    y = xy[2],
    size = pointSize,
    label = "id",
    customdata = "id"
  )

  # Only add colour and fill if we have valid variables
  if (!is.null(col_var) && col_var %in% colnames(df)) {
    aes_list$colour <- as.name(col_var)
  }
  if (!is.null(fill_var) && fill_var %in% colnames(df)) {
    aes_list$fill <- as.name(fill_var)
  }

  # Create plot
  p1 <- ggplot(df, aes_list) + geom

  if(!is.null(scales))
    p1 = p1 + scales
  if(!is.null(guides))
    p1 = p1 + guides
  if(!is.null(ylab))
    p1 = p1 +  ylab
    p1 = p1 +  theme_bw() + theme(aspect.ratio = 1,
                         panel.grid = element_blank(),
                         axis.text = element_text(color = "black"),
                         strip.background = element_rect(fill = "white"),
                         legend.key.height = unit(0.8, "lines"))
  p1
}
