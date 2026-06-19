
#' Plot Scatter for Bioconductor Experiment
#'
#' This function generates a scatter plot for a Bioconductor Experiment object.
#'
#' @param x A Bioconductor Experiment object.
#' @param chs A character vector specifying the variables to be plotted on the x and y axes.
#' @param gate Optional gate object.
#' @param color_by A character specifying the variable used for color coding points (optional).
#' @param facet_by A character specifying the variable used for faceting the plot (optional).
#' @param bins Number of bins for the 2D histogram (default is 100).
#' @param assay The assay to use for plotting (default is "exprs").
#' @param label Label for points on the plot ("target", "channel", or "both").
#' @param zeros Logical, whether to exclude rows with all zero values (default is FALSE).
#' @param k_pal Color palette for plotting.
#' @param rowNames Optional character vector to specify row names for the plot.
#'
#' @return A ggplot2 scatter plot.
#'
#' @examples
#' \dontrun{
#' # Requires a SingleCellExperiment with assay data.
#' plotScatterBJ(sce, chs = c("FSC.A", "SSC.A"))
#' }
#'
#' @export
plotScatterBJ <- function (x, chs, gate = NULL, color_by = NULL, facet_by = NULL, bins = 100,
                           assay = "exprs", label = c("target", "channel", "both"), 
                           zeros = FALSE, k_pal = NULL, rowNames = NULL)
{
  label <- match.arg(label)
  # args <- as.list(environment())
  # CATALYST:::.check_args_plotScatter(args)
  # rowNames takes some time to calculate, this optimizes for time
  if(is.null(rowNames)){
    m <- rownames(x)
  } else {
    m <- rowNames
  }
  if(is.null(chs) & !is.null(gate)){
    
  }
  # c <- channels(x)
  # browser()
  i = match(chs, m, nomatch = 0)
  if (all(i == 0)) i = NULL
  
  # i <- lapply(unique(m, c), function(u) {
  #   cat(file = stderr(), u, "\n")
  #   i <- match(chs, m, nomatch = 0)
  #   if (all(i == 0)) 
  #     NULL
  #   else i
  # })
  # i <- unlist(i) %>% unique()
  # y <- 
  y <- assay(x[i, , drop = FALSE], assay)
  chx <- CATALYST::channels(x)
  nms <- switch(label, target = m, channel = chx, both = ifelse(chx == m, chx, paste(chx, m, sep = "-")))
  chs[i != 0] <- rownames(y) <- nms[i]
  if (isTRUE(color_by %in% names(CATALYST::cluster_codes(x))))
    x[[color_by]] <- CATALYST::cluster_ids(x, color_by)
  cd <- cbind(SummarizedExperiment::colData(x), SingleCellExperiment::int_colData(x))
  df <- data.frame(t(as.matrix(y)), cd, check.names = FALSE, 
                   stringsAsFactors = FALSE)
  cd_vars <- intersect(names(cd), names(df))
  if (length(chs) > 2) {
    df <- reshape2::melt(df, id.vars = unique(c(chs[1], cd_vars)))
    facet <- "variable"
    ylab <- ylab(NULL)
    chs[2] <- "value"
  }  else {
    facet <- NULL
    ylab <- NULL
  }
  if (is.null(color_by)) {
    col_var <- guides <- NULL
    fill_var <- "..ncount.."
    if(ncol(x)<1000){
      scales <- scale_fill_gradientn(trans = "sqrt", colors = c("navy", "black"))
    } else{
      scales <- scale_fill_gradientn(trans = "sqrt", colors = c("navy", 
                                                                rev(brewer.pal(11, "Spectral"))))
    }
    geom <- geom_tile(stat = "bin2d", bins = bins, na.rm = TRUE, show.legend = FALSE)
  } else {
    fill_var <- NULL
    col_var <- sprintf("`%s`", color_by)
    geom <- geom_point(alpha = 0.2, size = 0.8, na.rm = TRUE)
    if (is.numeric(df[[color_by]])) {
      guides <- NULL
      scales <- scale_color_gradientn(colors = c("navy", 
                                                 rev(brewer.pal(11, "Spectral"))))
    }
    else {
      if (color_by %in% names(CATALYST::cluster_codes(x))) {
        if (is.null(k_pal)) k_pal <- CySA_default_cluster_cols()
        scales <- scale_color_manual(values = k_pal)
      }
      else scales <- NULL
      guides <- guides(col = guide_legend(override.aes = list(alpha = 1, 
                                                              size = 3)))
    }
  } 
  facet <- c(facet, facet_by)
  if (!is.null(facet)) {
    if (length(facet) == 1) {
      facet <- facet_wrap(facet)
    }
    else {
      facet <- facet_grid(cols = vars(!!rlang::sym(facet[1])),
                          rows = vars(!!rlang::sym(facet[2])))
    }
  }
  xy <- sprintf("`%s`", chs)
  if (!zeros) 
    df <- df[rowSums(df[, chs[c(1, 2)]] == 0) == 0, ]
  ggplot(df, aes(x = .data[[chs[1]]], y=.data[[chs[2]]], 
                 col = switch(is.null(color_by)+1,  .data[[color_by]], NULL)
                 # col = NULL,
                 # ,fill = {{fill_var}}
                 )) + 
    geom + scales + guides + facet + ylab + theme_bw() + 
    theme(aspect.ratio = 1, panel.grid = element_blank(), 
          axis.text = element_text(color = "black"), strip.background = element_rect(fill = "white"), 
          legend.key.height = unit(0.8, "lines"))
}
