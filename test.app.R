rootDir = "/pasteur/helix/projects/scBiomarkers/projects/julia/"

currentRoot = file.path(rootDir, "channel_rename", "keepChannels", "train", "indexed", "rmMargins", "noDoublets",  "cyEt", "transformed", "PeacoQC_results", "cyCombine", "scaleP.old", "sce.som.dim.40.rlen.200.mst.1.radius.40.seed.47.subsample400000")

cp = load(file = paste0(currentRoot,"/_data/sce.shiny.RData"))
devtools::load_all()

somCodes = metadata(sce)$SOM_codes

# raster plot of SOM
nDim = colnames(metadata(sce)$cluster_codes[,1,drop=F]) %>% stringr::str_sub(4) %>% as.integer() %>% sqrt()
somRasterData <- data.frame(x = rep(1:nDim, nDim),
                            y = rep(1:nDim, each = nDim)   # Create data frame for raster
)
somRasterData = cbind(somRasterData,somCodes)
colnames(somRasterData) = make.names(colnames(somRasterData))
somRasterObj <- raster::rasterFromXYZ(xyz = somRasterData)             # Convert data frame to raster object




# level1Markers = c("SSC.A","FSC.A","live.dead.A", "CD3.A", "CD16.A" , "CD19.A", "CD4.A", "CD14.A", "CD8.A")
dList = list(
  d1 = c("FSC.A", "SSC.A"),
  d2 = c("FSC.A", "FSC.H"),
  d3 = c("SSC.A", "live.dead.A"),
  d4 = c("CD3.A","SSC.A"),
  d5 = c("CD16.A", "CD19.A"),
  d6 = c("CD14.A", "CD14.A")
)

options(shiny.reactlog = TRUE)



devtools::load_all()

sh = clusterSelector(sce = sce, # main input has to contain:
                     sce_subsampled = sce_subsampled, # subsampled sce object
                     outputList = outputList, # list of named nodes
                     colTree = NULL, # Tree object to plot
                     dList = dList,
                     dend = dend,
                     dendTable = dendTable,
                     clusterPatientTable = clusterPatientTable,
                     somCodes = "SOM_codes", # SOM_codes.1
                     nPlots = 6,
                     somRasterData = somRasterData,
                     somRasterObj = somRasterObj)

outputList <- shiny::runApp(shinyApp(ui = sh[[1]], server = sh[[2]]))

outputList = get("outputList", envir = environment())

reactlog::reactlog_show()






