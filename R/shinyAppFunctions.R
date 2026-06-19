# functions used in shiny app

# highlight_df = function ----
highlight_df = function(x,y,rs, somCodesName = "SOM_codes", metaD = NULL){
  if (is.null(metaD)) stop("'metaD' must be provided")
  # cat(file = stderr(), paste(rs, collapse = ", "),"\n")
  data.frame(x=metaD[[somCodesName]][rs,x],
             y=metaD[[somCodesName]][rs,y],
             id=rs)
}



countBarPlotFunc <- function(rs, clusterPatientTable, cst, sce, outputList, groupsInput) {
  rs = as.integer(intersect(rs, colnames(clusterPatientTable)))
  # selected counts
  rSums = data.frame(counts = rowSums(clusterPatientTable[,rs,drop=F]))
  rSums$id = rownames(clusterPatientTable)

  # counts to compare to
  if(!cst=="none"){
    if(cst %in% colnames(S4Vectors::metadata(sce)$experiment_info)){
      # expInfo = S4Vectors::metadata(sce)$experiment_info
      numCols <- unlist(lapply(S4Vectors::metadata(sce)$experiment_info, is.numeric), use.names = FALSE)
      expInfo = S4Vectors::metadata(sce)$experiment_info[,numCols, drop=F]
      # eI = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]
      eI = apply(expInfo,2,as.numeric) %>% as.data.frame()
      rownames(eI) = S4Vectors::metadata(sce)$experiment_info[,"sample_id"]
      if(is.na(eI) %>% any()) {
        cat(file = stderr(), "some NAs produced 3\n")
        browser()
      }
      ctStats = eI[,cst]
      rSums = cbind(rSums, eI[rSums$id,cst])
      colnames(rSums) = c("counts", "id", cst)
      dfm <- reshape2::melt(rSums[,c('id',cst,'counts')],id.vars = 1)
      colnames(dfm) = c("id", "variable", "counts")
    }else if(cst %in% names(outputList)){
      rSums = cbind(rSums, rowSums(clusterPatientTable[,outputList[[cst]],drop=F]))
      colnames(rSums) = c("counts", "id", cst)
      dfm <- reshape2::melt(rSums[,c('id',cst,'counts')],id.vars = 1)
      colnames(dfm) = c("id", "variable", "counts")
    } else{
      cat(file = stderr(), "should not happen\n")
      browser()
    }
  } else {
    dfm <- rSums
    dfm$variable = "counts"
  }
  if("group" %in% colnames(S4Vectors::metadata(sce)$experiment_info)){
    dfm = merge(dfm, S4Vectors::metadata(sce)$experiment_info, by.x="id", by.y="sample_id")
    lvs=S4Vectors::metadata(sce)$experiment_info[order(S4Vectors::metadata(sce)$experiment_info$group),"sample_id"]
    # dfm = dfm %>% group_by(group)
    dfm$id = factor(dfm$id, levels=lvs)
  }
  dfm$groups = "other"
  if(!is_empty(groupsInput)){
    if(!is_empty(groupsInput$group1)){
      dfm$groups[dfm$id %in% groupsInput$group1] = "group1"
    }
    if(!is_empty(groupsInput$group2)){
      dfm$groups[dfm$id %in% groupsInput$group2] = "group2"
    }
  }
  dfm$groups = as.factor(dfm$groups)
  dfm = as_tibble(dfm)
  # save(file = "/pasteur/appa/scratch/bernd/Rtest2.Rdata",
  #      list = c("clusterPatientTable", "rs", "expInfo", "eI", "cst", "rSums", "dfm", "groupsInput"))
  # cp =load(file = "/pasteur/appa/scratch/bernd/Rtest2.Rdata")
  ggplot(dfm,aes(x = reorder(id,counts),y = counts)) +
    geom_bar(aes(fill = interaction(variable,groups)), stat = "identity", position = "dodge") +
    coord_flip() + guides(fill = guide_legend(title="groups"))
  # ggplot(dfm, aes(x=id, y=counts)) +
  #   geom_bar(stat = "identity") + coord_flip()
}


compute_relative_counts <- function(clusterPatientTable, rs, relativeToCol, expInfo, outputList) {
  switch(relativeToCol,
         "none" = rowSums(clusterPatientTable[, rs, drop = FALSE]) /
           rowSums(clusterPatientTable[, , drop = FALSE]) * 100,
         {
           if (relativeToCol %in% colnames(expInfo)) {
             rowSums(clusterPatientTable[, rs, drop = FALSE]) /
               expInfo[rownames(clusterPatientTable), relativeToCol] * 100
           } else if (relativeToCol %in% names(outputList)) {
             rowSums(clusterPatientTable[, rs, drop = FALSE]) /
               rowSums(clusterPatientTable[, outputList[[relativeToCol]], drop = FALSE]) * 100
           } else {
             stop("relativeToCol '", relativeToCol, "' not found")
           }
         }
  )
}

PercentBarPlotFunc <- function(sce, relativeToCol, clusterPatientTable, rs, outputList, group, groupsInput) {

  expInfo = metadata(sce)$experiment_info
  rownames(expInfo) = expInfo$sample_id
  expInfo = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]

  eI = apply(expInfo,2,as.numeric)
  # browser()
  # TODO make case
  # if(relativeToCol == "none"){
  #   rSums = data.frame(Percent = rowSums(clusterPatientTable[,rs,drop=F])/rowSums(clusterPatientTable[,,drop=F])*100)
  # } else {
  #   if (relativeToCol %in% colnames(metadata(sce)$experiment_info)){
  #     rSums = data.frame(Percent = rowSums(clusterPatientTable[,rs,drop=F])/expInfo[rownames(clusterPatientTable),relativeToCol]*100)
  #   } else {
  #     if(relativeToCol %in% names(outputList)){
  #       rSums = data.frame(Percent = rowSums(clusterPatientTable[,rs,drop=F])/rowSums(clusterPatientTable[,outputList[[relativeToCol]],drop=F])*100)
  #     } else{
  #       cat(file = stderr(), "ERROR\n")
  #     }
  #   }
  # }
  rSums <- data.frame(
    Percent = compute_relative_counts(
      clusterPatientTable, rs, relativeToCol, expInfo, outputList
    )
  )
  rSums$id = rownames(clusterPatientTable)

  if("group" %in% colnames(S4Vectors::metadata(sce)$experiment_info)){
    rSums = merge(rSums, S4Vectors::metadata(sce)$experiment_info, by.x="id", by.y="sample_id")
    rSums = rSums %>% group_by(group)
    rSums$id = factor(rSums$id, levels=rSums[order(rSums$group),"id"][[1]])
  }
  rSums$groups = "other"
  if(!is_empty(groupsInput)){
    if(!is_empty(groupsInput$group1)){
      rSums$groups[rSums$id %in% groupsInput$group1] = "group1"
    }
    if(!is_empty(groupsInput$group2)){
      rSums$groups[rSums$id %in% groupsInput$group2] = "group2"
    }
  }
  rSums$groups = as.factor(rSums$groups)
  rSums = as_tibble(rSums)
  # save(file = "/pasteur/appa/scratch/bernd/Rtest.Rdata",
  #      list = c("clusterPatientTable", "rs", "expInfo", "eI", "relativeToCol", "rSums", "groupsInput"))
  # cp = load(file = "/pasteur/appa/scratch/bernd/Rtest.Rdata")
  ggplot(rSums, aes(x=reorder(id,Percent), y=Percent, fill = groups)) +
    # geom_bar(stat = "identity") + coord_flip()
  geom_col() + coord_flip()

}



ggsomPlot <- function(pp1, plotIdx, rs, dimSelection, somCodesName = "SOM_codes", sce, metaD = S4Vectors::metadata(sce)){
  newData = highlight_df(dimSelection[[plotIdx]]$dims[1],dimSelection[[plotIdx]]$dims[2], rs, somCodesName, metaD = metaD)
  p3 = pp1 + geom_point(data=newData,
                        aes(x=`x`,y=`y`, customdata = rs),
                        color='red',
                        size=0.3)
  if (is.null(dimSelection[[plotIdx]]$xzoom[1])){
    # cat(file = stderr(), plotIdx, "xlim\n")

    p3 = p3 +
      xlim(c(dimSelection[[plotIdx]]$xlim[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$xlim[2] %>% as.numeric() )) +
      ylim(c(dimSelection[[plotIdx]]$ylim[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$ylim[2] %>% as.numeric() ))
  } else {
    # cat(file = stderr(), plotIdx, "zoomed\n")
    p3 = p3 +
      xlim(c(dimSelection[[plotIdx]]$xzoom[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$xzoom[2] %>% as.numeric() )) +
      ylim(c(dimSelection[[plotIdx]]$yzoom[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$yzoom[2] %>% as.numeric() ))
  }
  return(p3)
}

somPlot <- function(pp1, plotIdx, rs, colorbyGroups, showGroups, dimSelection = NULL, somCodesName = "SOM_codes", sce, metaD = S4Vectors::metadata(sce), outputList = list(), projectionDf = NULL){
  if(is.null(pp1)) return(NULL)

  # cp =load(file = "/pasteur/appa/scratch/bernd/dev.RData")
  # if(plotIdx==1){
  #   # browser()
  #   save(file = "/pasteur/appa/scratch/bernd/dev.RData", list = c("pp1", "plotIdx", "rs", "dimSelection", "colorbyGroups",
  #                                                                 "showGroups"))
  # }
  if(showGroups) {
    req(projectionDf)
    p3 = drawProjection(projectionDf, rs, colorbyGroups = colorbyGroups, sce = sce, outputList = outputList)
  } else {
    p3 = pp1 + geom_point(data=highlight_df(dimSelection[[plotIdx]]$dims[1],dimSelection[[plotIdx]]$dims[2], rs, somCodesName, metaD = metaD),
                          aes(x=x,y=y, customdata=rs),
                          color='red',
                          size=0.3)
  }

  if(is.null(p3)) return(NULL)
  if (is.null(dimSelection[[plotIdx]]$xzoom[1])){
    p3 = p3 +
      xlim(c(dimSelection[[plotIdx]]$xlim[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$xlim[2] %>% as.numeric() )) +
      ylim(c(dimSelection[[plotIdx]]$ylim[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$ylim[2] %>% as.numeric() ))
  } else {
    # cat(file = stderr(), plotIdx, "zoomed\n")
    p3 = p3 +
      xlim(c(dimSelection[[plotIdx]]$xzoom[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$xzoom[2] %>% as.numeric() )) +
      ylim(c(dimSelection[[plotIdx]]$yzoom[1] %>% as.numeric() ,
             dimSelection[[plotIdx]]$yzoom[2] %>% as.numeric() ))
  }
  # browser()
  if(plotIdx>5) plotIdx = 6


  ggplotly(p3, source = paste0("somData", plotIdx), tooltip = "") %>%
    layout(showlegend = F, dragmode = "select") %>%
    config(renderer = "webgl") %>%
    event_register("plotly_selected") %>%
    event_register("plotly_relayout")
}




tsneFunc <- function(dimRedSelection, perplexity, sce, somCodesName = "SOM_codes") {
  seed = 1
  dimRedCols = dimRedSelection
  set.seed(seed)
  Rtsne::Rtsne(S4Vectors::metadata(sce)[[somCodesName]][,dimRedCols], perplexity = perplexity)
}



plotViolinFunc <- function(sce, somCodesName = "SOM_codes", upsetSelection, outputList, violinSelection) {
  markers = colnames(sce@metadata[[somCodesName]])
  if(length(upsetSelection)<3)
    upsetSelection = names(outputList)
  markers = intersect(markers, violinSelection)
  data = data.frame(
    somNode = factor(levels = upsetSelection),
    marker = factor(levels = markers),
    expr = double(),
    grpName = factor(levels = upsetSelection)

  )
  for( na in upsetSelection){
    wide = sce@metadata[[somCodesName]][outputList[[na]],markers,drop=FALSE] %>% as.data.frame()
    if(nrow(wide)<2)next()
    if(any(outputList[[na]] ==0)){
      outputList[[na]] = outputList[[na]][-which(outputList[[na]] ==0)]
    }

    wide$somNode = factor(outputList[[na]], levels = 1:nrow(sce@metadata[[somCodesName]]))
    long = gather(wide, marker,expr,-somNode, factor_key=TRUE)
    long$grpName = na
    data = rbind(data, long)
  }
  if(nrow(data)>10){
    nb.cols <- length(unique(data$marker))
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

    p = ggplot(data, aes(factor(marker), expr, fill = marker)) +
      geom_violin(scale = "width", adjust = 1, trim = TRUE) +
      scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
        c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
      facet_grid(rows = vars(grpName), scales = "free", switch = "y") +
      theme_cowplot(font_size = 12) +
      theme(legend.position = "none", panel.spacing = unit(0, "lines"),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = NA, color = "black"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.y.left = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = mycolors) +
      ggtitle("Marker on x-axis") + xlab("Marker") + ylab("Expression Level")
    return(p)
  } else{
    message("please check that all upsetSelection are in outputList: no violin data to plot")
    return(NULL)
  }
}



# Build a data frame of projected SOM coordinates plus SOM_stats, cached per view.
buildProjectionDf <- function(pp1, plotIdx, dimSelection, sce) {
  req(pp1, dimSelection)
  pg <- ggplot_build(pp1)
  df <- pg$data[[1]][, c("x", "y", "label")]
  df2 <- dplyr::left_join(S4Vectors::metadata(sce)$SOM_stats, df,
                          by = dplyr::join_by(id == label))
  df2 <- df2[order(df2$id), ]
  df2[is.na(df2)] <- 0
  ch_names <- dimSelection[[plotIdx]]$dims
  colnames(df2)[colnames(df2) == "x"] <- ch_names[1]
  colnames(df2)[colnames(df2) == "y"] <- ch_names[2]
  # Coordinate columns must be first so drawProjection can locate them
  coord_idx <- match(ch_names, names(df2))
  other_idx <- setdiff(seq_along(df2), coord_idx)
  df2 <- df2[, c(coord_idx, other_idx), drop = FALSE]
  df2
}

drawProjection <- function(df, rs, colorbyGroups, sce, outputList = list()){
  colN = names(df)[1:2]
  df$cluster = df$id
  # N, mean, thrdQu, max are already in df from SOM_stats
  nGrps = 1
  if(length(colorbyGroups)<1){
    colGrp = "lightblue"
    sl = FALSE
  } else {
    df$colGrp = ""
    for(cg in colorbyGroups){
      df$colGrp[df$cluster %in% outputList[[cg]]] = paste(df$colGrp[df$cluster %in% outputList[[cg]]],cg)
    }
    df$colGrp[df$colGrp == ""] = "other"
    df$colGrp = factor(df$colGrp)
    sl = TRUE
    nGrps = length(levels(df$colGrp))
    nb.cols <- nGrps + 1
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
  }
  # browser()
  p3 = ggplot(data = df,
              aes(x=.data[[colN[1]]], y=.data[[colN[2]]],
                  text = paste("cluster: ",cluster,"<br>N cells: ",N,"<br>mean: ",
                               format(mean,digits=2),"<br>3rd Q: ",
                               format(thrdQu,digits=2),"<br>max: ",
                               format(max,digits=2)),
                  customdata = seq(nrow(df))
              )
  )
  if(nGrps==1){
    p3 = p3 +geom_point(show.legend = sl, color="lightblue")
  } else{
    p3 = p3 + geom_point(show.legend = sl, aes(color=colGrp)) +
      scale_color_manual(values = mycolors[seq_len(nGrps)])
  }

  p3 = p3 +
    geom_point(data=df[rs,],
               aes(x=.data[[colN[1]]],
                   y=.data[[colN[2]]],
                   customdata=rs),
               color='red',
               size=0.3) +
    theme_minimal() +
    theme(legend.position="bottom") +
    guides(colour = guide_legend(nrow = as.integer(nGrps/3+1),title = ""))

  return(p3)
}


plotViolin2Func <- function(sce, somCodesName = "SOM_codes", violinSelection, upsetSelection, outputList) {

  markers = colnames(sce@metadata[[somCodesName]])
  markers = intersect(markers, violinSelection)
  if(length(upsetSelection)<3)
    upsetSelection = names(outputList)
  data = data.frame(
    somNode = factor(levels = upsetSelection),
    marker = factor(levels = markers),
    expr = double(),
    grpName = factor(levels = upsetSelection)

  )
  for( na in upsetSelection){
    if(length(outputList[[na]])==0)next()
    wide = sce@metadata[[somCodesName]][outputList[[na]],markers,drop=FALSE] %>% as.data.frame()
    if(nrow(wide)<2)next()
    if(any(outputList[[na]] ==0)){
      outputList[[na]] = outputList[[na]][-which(outputList[[na]] ==0)]
    }
    wide$somNode = factor(outputList[[na]], levels = 1:nrow(sce@metadata[[somCodesName]]))
    long = gather(wide, marker,expr,-somNode, factor_key=TRUE)
    long$grpName = na
    data = rbind(data, long)
  }
  nb.cols <- length(unique(data$grpName))
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

  p = ggplot(data, aes(factor(grpName), expr, fill = grpName)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(marker), scales = "free", switch = "y") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = mycolors) +
    ggtitle("grp name on x-axis") + xlab("Groups") + ylab("Expression Level")
  p

}



upsetPlotFunc <- function(upsetSelection, outputList, sce) {

  if(length(upsetSelection)>31) upsetSelection = upsetSelection[1:31]
  cm = ComplexHeatmap::make_comb_mat(outputList[upsetSelection])
  ncells = rep(0,length(ComplexHeatmap::comb_name(cm)))
  grpCells = rep(0,length(upsetSelection))
  for(cidx in seq(ComplexHeatmap::comb_name(cm))){

    ncells[cidx] = S4Vectors::metadata(sce)$SOM_stats[as.integer(ComplexHeatmap::extract_comb(cm, ComplexHeatmap::comb_name(cm)[cidx])),"n"] %>% sum()
  }
  for(olIdx in seq(upsetSelection)){
    grpCells[olIdx] =  S4Vectors::metadata(sce)$SOM_stats[as.integer(outputList[[upsetSelection[olIdx]]]),"n"] %>% sum()
  }
  names(ncells) = ComplexHeatmap::comb_name(cm)
  names(grpCells) = upsetSelection
  grpSoms = lapply(outputList[upsetSelection],length) %>% unlist()
  top_ha = ComplexHeatmap::HeatmapAnnotation(
    "cell #" = ComplexHeatmap::anno_barplot(ncells, add_numbers = TRUE,
                            gp = grid::gpar(fill = "black"), width = grid::unit(4, "cm")),
    "som #" = ComplexHeatmap::anno_barplot(ComplexHeatmap::comb_size(cm), add_numbers = TRUE,
                           gp = grid::gpar(fill = "black"), width = grid::unit(4, "cm")),
    gap = grid::unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0
  )

  side_ha = ComplexHeatmap::rowAnnotation(
    "som #" = ComplexHeatmap::anno_barplot(grpSoms, add_numbers = TRUE,
                           gp = grid::gpar(fill = NULL), width = grid::unit(3, "cm")),

    "cell #" = ComplexHeatmap::anno_barplot(grpCells, add_numbers = TRUE,
                            gp = grid::gpar(fill = NULL), width = grid::unit(3, "cm")),
    gap = grid::unit(2, "mm"),  annotation_name_rot = 0
  )
  ComplexHeatmap::UpSet(cm, comb_order = order(ComplexHeatmap::comb_degree(cm), -ComplexHeatmap::comb_size(cm)),
        top_annotation = top_ha, right_annotation = side_ha)

}

