utils::globalVariables(c(
  ".", ".data", "..ncount..", "counts", "groups", "id", "marker", "expr",
  "grpName", "Percent", "n", "value", "variable", "cluster", "colGrp",
  "x", "y", "customdata", "label", "tsne1", "tsne2", "umap1", "umap2",
  "pc1", "pc2", "rowElement", "cellCounts", "cellPercentages",
  "dend", "dendTable", "clusterPatientTable", "somRasterData",
  "somRasterObj", "sce_subsampled", "metaD", "sceRN", "sceCN",
  "chx", "nNsub", "sce_subsampledRN", "df", "dfPlot",
  "outputList", "dimSelection", "activePlot", "rsUsed",
  "inputClusterNumber", "violinPlotSelection", "groupsInput",
  "countBarPlot", "PercentBarPlot", "dendPlot", "tsnePlot",
  "umapPlot", "pcaPlot", "scatterPlot", "somRasterPlot",
  "vlnPlot", "VlnPlot2", "upSetPlot", "selectedUpdate", "selectedUpdate2",
  "sample2PlotDb", "choicesRV",
  "sample_id", "group", "val", "somNode", "N", "thrdQu"
))

#' Cluster Selector Shiny Application
#'
#' Creates an interactive Shiny dashboard for selecting and visualizing
#' clusters from a SingleCellExperiment object.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing the full dataset.
#' @param sce_subsampled A subsampled \code{SingleCellExperiment} for
#'   performance-sensitive plots.
#' @param outputList A named list of cluster groupings.
#' @param colTree Optional collapsible tree object.
#' @param dList List of marker pairs for 2D plots.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param nPlots Number of 2D SOM plots to display.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somRasterObj Raster object for SOM visualization.
#' @param env Environment used to store mutable state (legacy argument).
#'
#' @return A list with two elements: \code{ui} and \code{server}.
#'
#' @export
clusterSelector <- function(sce, # main input has to contain:
                            sce_subsampled, # subsampled sce object
                            outputList=list(), # list of named nodes
                            colTree=NULL, # Tree object to plot
                            dList,
                            dend,
                            dendTable,
                            clusterPatientTable,
                            somCodesName = "SOM_codes", # SOM_codes.1
                            nPlots = 6,
                            somRasterData,
                            somRasterObj,

                            env = environment()) {

  for(idx in seq(dList)){
    assign(paste0("d",idx,".1"), dList[[idx]][1])
    assign(paste0("d",idx,".2"), dList[[idx]][2])
  }

  if(!is.table(clusterPatientTable)){
    stop("clusterPatientTable not a table")
  }

  if(!"cluster_id" %in% names(dimnames(clusterPatientTable))){
    stop("cluster_id not in colnames(clusterPatientTable)")
  }

  rnSCE = rownames(sce)

  jscode <- "shinyjs.closeWindow = function() { window.close(); }"

  if("outputList" %in% ls(envir = env)){
    outputList = get("outputList", envir = env)
  }
  outputList[["Rest"]] = c()
  outputList[["Rest"]] = as.integer(levels(sce$cluster_id)[!levels(sce$cluster_id) %in% unique(unlist(outputList))])
  # outputList = append(outputList, list(cList ))
  # names(outputList)[length(outputList)] = cName
  for(na in names(outputList)){
    if(length(outputList[[na]])==0) outputList[[na]] = NULL
  }
  metaD = S4Vectors::metadata(sce)
  colsUsed = metaD$map$colsUsed
  # only necessary data for raster plots.
  metaD$map$colsUsed
  # all SOM cols have to be included
  missing_cols <- setdiff(metaD$map$colsUsed, colnames(somRasterData))
  if (length(missing_cols) > 0) {
    stop("somRasterData is missing SOM columns: ", paste(missing_cols, collapse = ", "))
  }
  somRasterData = somRasterData[,c("x", "y", metaD$map$colsUsed)]
  somRasterObj = somRasterObj[[metaD$map$colsUsed]]

  # Pre-compute per-channel axis limits once from the full assay matrix.
  # This avoids scanning assays(sce)[[1]] every time dimSelection changes.
  assay_mat <- SummarizedExperiment::assays(sce)[[1]]
  channelLimits <- lapply(rownames(sce), function(ch) {
    vals <- assay_mat[ch, ]
    c(min = min(vals), max = max(vals))
  })
  names(channelLimits) <- rownames(sce)

  assign(x = "outputList", value = outputList, envir = env)

  # nPlots = 15
  # nPlots = 6

  {
    ## UI ----
    ui <- dashboardPage(
      dashboardHeader(title = "Cluster Selector"),

      sidebar = dashboardSidebar(
        shinyjs::useShinyjs(),
        shinyjs::extendShinyjs(text = jscode, functions = c("closeWindow")),

        sidebarMenu(id="tabs",
                    menuItem("Parameters to adjust", tabName="parameters", icon=icon("line-chart")
                    )
        ),
        hr(),
        ### parameters ----
        conditionalPanel("input.tabs == 'parameters'",
                         fluidRow(
                           column(1),
                           column(10,
                                  htmltools::tags$p("The node you most recently clicked:"),
                                  verbatimTextOutput("str")),
                           # selectInput("d1.1",
                           #             "select X (plot 1)",
                           #             choices = rownames(sce),
                           #             multiple = F,selectize = T),

                           shinydashboardPlus::box(
                             title = "dim modifications", solidHeader = TRUE, width = 12,
                             collapsible = TRUE, collapsed = FALSE,
                             uiOutput("dimUI")
                           ),
                           # uiOutput("axesUI"),
                           {
                             choices = lapply(dList, FUN = function(x) paste(x, collapse = " - ")) %>% unlist() %>% unname()

                             selectizeInput(
                               inputId = "d2Axes",
                               label = "Choose a plot index :",
                               choices = choices,
                               selected = choices[1:6],
                               multiple = TRUE,
                               size = 6
                             )
                           },
                           actionButton("applyDimSelection", "apply dim selection")
                           ,
                           # radioButtons("d2Axes",
                           #              "Choose a plot index :",
                           #              choices = lapply(dList, FUN = function(x) paste(x, collapse = " - ")) %>% unlist() %>% unname(),
                           #              selected = 1,
                           #              inline = FALSE
                           # ),
                           radioButtons(
                             "selectMode",
                             "select mode",
                             choices = c("view", "remove others", "add", "remove"),
                             selected = "view",
                             inline = FALSE
                           ),
                           selectInput("samples2plot",
                                       paste0("samples to plot"),
                                       choices = levels(sce$sample_id),
                                       selected = levels(sce$sample_id),
                                       multiple = T,selectize = T),
                           textInput("clusterNumbers",
                                     "cluster numbers",
                                     value = "1",
                                     width = NULL,
                                     placeholder = NULL),
                           actionButton("applyclusterNumbers", "apply cluster numbers"),
                           textInput("clusterName",
                                     "name selection",
                                     value = "",
                                     width = NULL,
                                     placeholder = NULL),
                           actionButton("applyName", "apply name"),
                           selectInput("groupRM",
                                       "select to remove",
                                       choices = names(outputList),
                                       multiple = T,selectize = T),
                           actionButton("rmGroups", "remove groups"),
                           selectInput("clusterNameSelect",
                                       "select named",
                                       choices = names(outputList),
                                       multiple = T,selectize = T),
                           # checkboxInput("showPoints","show individual points", FALSE),
                           selectInput("clusterNameRM",
                                       "select to remove",
                                       choices = names(outputList),
                                       multiple = F,selectize = T),
                           actionButton("rmGrp", "remove"),
                           downloadButton('downloadPlots', 'Download Plots'),
                           actionButton("close", "Close window")
                         )
        )
      ),
      # Sidebar with a select input for the root node
      ### Body ----

      body = dashboardBody(
        ### First row ----
        fluidRow(
          if(!is.null("colTree")){
            column(width = 3,
                   shinydashboardPlus::box(
                     title = "interactive tree", solidHeader = TRUE, width = 12, status = "primary",
                     collapsible = TRUE, collapsed = TRUE,
                     collapsibleTree::collapsibleTreeOutput("plot"))
            )
          },
          column(width = 6,
                 shinydashboardPlus::box(
                   title = "2D plot", solidHeader = TRUE, width = 12, status = "primary",
                   collapsible = TRUE, collapsed = TRUE,
                   plotly::plotlyOutput("scatter") %>% shinyjqui::jqui_resizable()
                 )
          ),
          column(width = 3,
                 shinydashboardPlus::box(
                   title = "dendrogram", solidHeader = TRUE, width = 12, status = "primary",
                   collapsible = TRUE, collapsed = TRUE,

                   shiny::plotOutput("dend") %>% shinyjqui::jqui_resizable()
                 )
          )
        ),
        ### som Data ----
        shinydashboardPlus::box(
          title = "som 2D plots", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = TRUE, collapsed = TRUE,
          fluidRow(column(width = 2,
                          checkboxInput("showGroups", "color groups",value = FALSE))),
          fluidRow(column(width = 4,
                          plotly::plotlyOutput("somData1") %>% shinyjqui::jqui_resizable() ),
                   column(width = 4,
                          plotly::plotlyOutput("somData2") %>% shinyjqui::jqui_resizable()),
                   column(width = 4,
                          plotly::plotlyOutput("somData3") %>% shinyjqui::jqui_resizable())
          ),
          if(nPlots>3){
            fluidRow(
              column(width = 4,
                     plotly::plotlyOutput("somData4") %>% shinyjqui::jqui_resizable()),
              if(nPlots>4){
                column(width = 4,
                       plotly::plotlyOutput("somData5") %>% shinyjqui::jqui_resizable())
              },
              if(nPlots>5){
                column(width = 4,
                       plotly::plotlyOutput("somData6") %>% shinyjqui::jqui_resizable()
                )}

            )
          },
          fluidRow(column(width = 8,
                          selectInput("colorbyGroups",
                                      paste0("Select groups to color by"),
                                      choices = names(outputList),
                                      multiple = T,selectize = T),
                          {
                            cn = setdiff(colsUsed, c("label", "clusterid"))
                            selectInput("dimRedSelection",
                                        paste0("Select markers to use for dim. Red."),
                                        choices = cn,
                                        selected = cn,
                                        multiple = T,selectize = T)
                          }
          ),
          column(width = 2,
                 numericInput("perplexity",
                              "Perplexity",
                              value = 30, min = 1, max = 500),
                 checkboxInput("showlegend","show legend",value = F)
          ),
          column(width = 2,
                 numericInput("n_neighbors",
                              "n_neighbors",
                              value = 4, min = 2, max = 500)
          )

          ),
          fluidRow(column(width = 4,
                          plotly::plotlyOutput("tsne") %>% shinyjqui::jqui_resizable()),
                   column(width = 4,
                          plotly::plotlyOutput("umap") %>% shinyjqui::jqui_resizable()),
                   column(width = 4,
                          plotly::plotlyOutput("pca") %>% shinyjqui::jqui_resizable())
          )),
        ### Stats ----
        shinydashboardPlus::box(
          title = "Stats", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = TRUE, collapsed = TRUE,
          fluidRow(column(width = 4,
                          verbatimTextOutput("somClusters")),
                   column(width=8,
                          fluidRow(column(width=4,
                                          selectInput("compareStatsTo",
                                                      "select comparison Stats",
                                                      choices = {
                                                        numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
                                                        expInfo = metaD$experiment_info[,numCols, drop=F]
                                                        # expInfo = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]
                                                        eI = apply(expInfo,2,as.numeric)
                                                        if(is.na(eI) %>% any()) {
                                                          cat(file = stderr(), "some NAs produced\n")
                                                          browser()
                                                        }
                                                        c("none",colnames(eI))
                                                      },
                                                      selected = "none",
                                                      multiple = F,selectize = T)),
                                   column(width=4,
                                          selectInput("relativeTo",
                                                      "Parent population for percentages",
                                                      choices = {
                                                        numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
                                                        expInfo = metaD$experiment_info[,numCols, drop=F]
                                                        eI = apply(expInfo,2,as.numeric)
                                                        if(is.na(eI) %>% any()) {
                                                          cat(file = stderr(), "some NAs produced")
                                                          browser()
                                                        }
                                                        c("none",colnames(eI))
                                                      },
                                                      selected = "none",
                                                      multiple = F,selectize = T)),
                                   column(width=4,
                                          numericInput("singleNode",
                                                       "single SOM node counts",
                                                       min = 1,
                                                       max = nrow(metaD[[somCodesName]]),
                                                       value = 1, step =1))
                          ),
                          fluidRow(
                            column(width=4,
                                   selectInput("groupsVar",
                                               "Select groups variable for t-test",
                                               choices = {
                                                 factCols <- unlist(lapply(metaD$experiment_info, is.factor), use.names = FALSE)
                                                 # expInfo = metaD$experiment_info[,factCols, drop=F]
                                                 # expInfo = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]
                                                 # eI = apply(expInfo,2,as.numeric)
                                                 # if(is.na(eI) %>% any()) {
                                                 #   cat(file = stderr(), "some NAs produced\n")
                                                 #   browser()
                                                 # }
                                                 c("none",colnames(metaD$experiment_info)[factCols])
                                               },
                                               selected = "none",
                                               multiple = F,selectize = T)
                            ),
                            column(width=4,
                                   selectInput("group1",
                                               "group 1",
                                               choices = c(),
                                               selected = "",
                                               multiple = T, selectize = T
                                   )
                            ),
                            column(width=4,
                                   selectInput("group2",
                                               "group 2",
                                               choices = c(),
                                               selected = "",
                                               multiple = T, selectize = T
                                   )
                            )
                          ),
                          fluidRow(
                            shinydashboardPlus::box(title = "Stats", solidHeader = TRUE, width = 12, status = "primary",
                                                    collapsible = TRUE, collapsed = TRUE,
                                                    fluidRow(
                                                      DT::DTOutput("cellCounts") %>% shinyjqui::jqui_resizable(),
                                                      shiny::verbatimTextOutput("cellPercentages") %>% shinyjqui::jqui_resizable(),
                                                    ),
                                                    fluidRow(column(width = 6,

                                                                    shiny::plotOutput("CountBar") %>% shinyjqui::jqui_resizable()),
                                                             column(width = 6,

                                                                    shiny::plotOutput("PercentBar") %>% shinyjqui::jqui_resizable()
                                                             )
                                                    )
                            )),
                          fluidRow(
                            shinydashboardPlus::box(title = "t-Test", solidHeader = TRUE, width = 12, status = "primary",
                                                    collapsible = TRUE, collapsed = TRUE,

                                                    fluidRow(
                                                      column(width = 12,
                                                             verbatimTextOutput("ttestResult")
                                                      )
                                                    )
                            )
                          )
                   )
          )),
        ### som Plots ----
        shinydashboardPlus::box(
          title = "som plots", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = TRUE, collapsed = TRUE,
          fluidRow(column(width = 12,
                          shiny::plotOutput("somRaster", height = "1200px") %>% shinyjqui::jqui_resizable())),
          fluidRow(column(width = 3,
                          plotly::plotlyOutput("somRasterSelect", height = "400px") %>% shinyjqui::jqui_resizable()))
        ),
        ### Violin plots  ----
        shinydashboardPlus::box(
          id = "violinBox",
          title = "violin plots", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = TRUE, collapsed = TRUE,
          fluidRow(column(width = 12),
                    selectInput("violinSelection",
                                paste0("Select markers to show"),
                                choices = colnames(sce@metadata[[somCodesName]]),
                                selected = colsUsed,
                                multiple = T,selectize = T)),
          fluidRow(column(width = 6,
                          shiny::plotOutput("VlnPlot") %>% shinyjqui::jqui_resizable()),
                   column(width = 6,
                          shiny::plotOutput("VlnPlot2") %>% shinyjqui::jqui_resizable())
          )
        ),
        ### Upset Plots ----
        shinydashboardPlus::box(
          id = "upsetBox",
          title = "UpSet plot", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = TRUE, collapsed = TRUE,
          fluidRow(column(width = 12),
                   selectInput("upsetSelection",
                               paste0("Select groups to show"),
                               choices = names(outputList),
                               selected = names(outputList),
                               multiple = T,selectize = T)),
          fluidRow(column(width = 12,
                          shiny::plotOutput("UpSet") %>% shinyjqui::jqui_resizable())
          )
        )
        # ,
        # fluidRow(column(width = 4,
        #                 plotlyOutput("scatterPoints"))
        #          # ,
        #          # column(width = 4,
        #          #        )
        # )
      )

      # Show a tree diagram with the selected root node


    )

    #  server  ----
    server <- function(input, output, session) {
      rsUsed <- reactiveVal(c(1))
      triggerRedraw <- reactiveVal(1)
      selectedPoints = reactiveVal(NULL)
      activePlot <- reactiveVal(1)
      df = as.data.frame(t(assays(sce_subsampled)[[1]]))
      colnames(df) <- make.names(rownames(sce_subsampled))
      dfPlot = df
      # dfall = data.frame(t(assays(sce)[[1]]))
      # rN = rownames(sce)
      metaD = S4Vectors::metadata(sce)
      rN = colnames(metaD[[somCodesName]])
      nNsub = colnames(S4Vectors::metadata(sce_subsampled)[[somCodesName]])
      rownames(sce_subsampled)
      sce_subsampledRN = rownames(sce_subsampled)
      sceRN = rownames(sce)
      sceCN = colnames(sce)
      chx = channels(sce)
      if(!is.null("colTree")){
        output$plot <- collapsibleTree::renderCollapsibleTree({
          colTree
        })
      }

      # force redraw if selected group is used
      selectedUpdate2 <- reactiveVal(value=0)

      # reactiveValues replacement for external env$outputList
      rv <- reactiveValues(outputList = outputList)

      # Keep external env in sync for callers that still read env$outputList
      shiny::observe({
        assign(x = "outputList", value = rv$outputList, envir = env)
      })

      shiny::observe({
        rs = rsUsed_d()
        if("selected" %in% input$colorbyGroups){
          isolate(selectedUpdate2(selectedUpdate2() +1))
        }
      })

      # dim UI ----
      #

      choicesRV <- reactiveValues(trigger=1)

      observe({
        choicesRV$trigger
        cat(file=stderr(), "choicesRV triggered\n")
      })

      output$dimUI <- renderUI({
        cat(file = stderr(), "dimUI\n")
        cat(file = stderr(), "dimUI trigger: ", choicesRV$trigger, "\n")
        triggered = choicesRV$trigger
        if(!is.null(triggered)){
          cat(file = stderr(), "dimUI trigger: ", choicesRV$trigger, "\n")
        }
        tagsX = list()
        if(!base::exists("d1.1")){ # lets assume that if one is set all are set
          ridx = 1
          for(idx in 1:nPlots){
            assign(paste0("d",idx,".1"), rownames(sce)[ridx])
            ridx = ridx + 1
            assign(paste0("d",idx,".2"), rownames(sce)[ridx])
            ridx = ridx + 1
          }
        }
        for(idx in 1:nPlots){
          tagsX[[length(tagsX)+1]] = fluidRow(column(width = 6,
                                                     selectInput(paste0("d",idx,".1"),
                                                                 paste0("select X (plot ", idx, ")"),
                                                                 choices = rownames(sce),
                                                                 selected = get(paste0("d",idx,".1")),
                                                                 multiple = F,selectize = T),
          ),
          column(width = 6,
                 selectInput(inputId = paste0("d",idx,".2"),
                             label = paste0("select Y (plot ", idx, ")"),
                             choices = rownames(sce),
                             selected = get(paste0("d",idx,".2")),
                             multiple = F,selectize = T),
          )
          )

        }

        tagList(tagsX)
      })


      shiny::observeEvent(input$applyDimSelection, {
        cat(file = stderr(), "============inputList observer\n")

        inputList = input$d2Axes
        newVals = lapply(inputList[1:6], FUN = function(x) str_split(string = x,pattern = " - ")) %>% unlist()
        ridx = 1
        for(idx in 1:nPlots){
          assign(paste0("d",idx,".1"), newVals[ridx])
          ridx = ridx + 1
          assign(paste0("d",idx,".2"), newVals[ridx])
          ridx = ridx + 1
        }
        for(idx in 1:nPlots){

          updateSelectInput(session = session,
                            inputId = paste0("d",idx,".1"),
                            selected = get(paste0("d",idx,".1")),
          )
          updateSelectInput(session = session,
                            inputId = paste0("d",idx,".2"),
                            selected = get(paste0("d",idx,".2")),
          )


        }
        cat(file = stderr(), "triggering choicesRV\n")
        choicesRV$trigger = choicesRV$trigger + 1
      })


# update Radio buttons d2Axes1 ----
# lapply(1:nPlots, function(x){
#   shiny::observeEvent(input[[paste0("d",x,".2")]],{
#     choices = list()
#     for(idx in 1:nPlots){
#       if(input[[paste0("d", idx, ".1")]] %in% rnSCE ){
#         # choices[[length(choices) + 1]] = paste0(input[[paste0("d", idx, ".1")]], "/", input[[paste0("d", idx, ".2")]])
#         choices[[paste0(input[[paste0("d", idx, ".1")]], "/", input[[paste0("d", idx, ".2")]])]] = idx
#       }
#
#     }
#
#     cat(file = stderr(), "here i am", x,"\n")
#     oldVal = input$d2Axes1
#     updateRadioButtons(inputId = "d2Axes1", choices = choices, selected = oldVal)
#   })
# })

# UI Choose a plot index ----
output$axesUI <- renderUI({
  cat(file = stderr(), "UI\n")
  if(is.null(input[[paste0("d", idx, ".1")]])) return(NULL)
  choices = list()
  for(idx in 1:nPlots){
    if(input[[paste0("d", idx, ".1")]] %in% rnSCE ){
      choices[[paste0(input[[paste0("d", idx, ".1")]], "/", input[[paste0("d", idx, ".2")]])]] = idx
    }

  }
  # browser()
  radioButtons("d2Axes1",
               "Choose a plot index :",
               choices = choices,
               selected = 1,
               inline = FALSE
  )
}
)

# observe clusterNameSelect ----
# print clusters based on selection
shiny::observeEvent(input$clusterNameSelect,{
  outputList = rv$outputList
  listNames = names(outputList) %in% input$clusterNameSelect
  combinedSoms =  outputList[listNames] %>% unlist() %>% unique()
  updateTextInput(session, "clusterNumbers", value  = paste(combinedSoms, collapse=", "))
})

# observeEvent groupsVar ----
shiny::observeEvent(input$groupsVar,{
  # browser()
  groupsVar = input$groupsVar
  if(! input$groupsVar %in% colnames(metaD$experiment_info)) return(NULL)
  levs = levels(metaD$experiment_info[,input$groupsVar])
  updateSelectInput(session = session, inputId = "group1", choices = levs)
  updateSelectInput(session = session, inputId = "group2", choices = levs)
})

# observe group ----
shiny::observe({
  if(! input$groupsVar %in% colnames(metaD$experiment_info)) return(NULL)
  levs = levels(metaD$experiment_info[,input$groupsVar])
  grp1 = input$group1
  grp2 = isolate(input$group2)
  levs = levs[!levs %in% grp1]
  updateSelectInput(session = session, inputId = "group2", choices = levs, selected = grp2)
  # save(file = "dev2.RData", list = c("groupsVar", "levs", "grp1"))
  #
})
shiny::observe({
  if(! input$groupsVar %in% colnames(metaD$experiment_info)) return(NULL)
  levs = levels(metaD$experiment_info[,input$groupsVar])
  grp2 = input$group2
  grp1 = isolate(input$group1)
  levs = levs[!levs %in% grp2]
  updateSelectInput(session = session, inputId = "group1", choices = levs, selected = grp1)
  # save(file = "dev2.RData", list = c("groupsVar", "levs", "grp1"))
  #
})
# delayed group input for t-test ----
# returns sample_ids from potentially other metadata factorials
groupsInput <- reactive({
  if(! input$groupsVar %in% colnames(metaD$experiment_info)) return(NULL)
  cat(file = stderr(), "delayed Groups input: \n")
  grp1 = metaD$experiment_info[
    metaD$experiment_info[,input$groupsVar] %in% input$group1,"sample_id"]
  grp2 = metaD$experiment_info[
    metaD$experiment_info[,input$groupsVar] %in% input$group2,"sample_id"]
  list(group1 = grp1, group2 = grp2 )
}) %>% debounce(1000)

# t.test ----
shiny::observe({
  cat(file = stderr(), "t.test\n")
  empty = FALSE
  groupsVar = input$groupsVar
  if(! input$groupsVar %in% colnames(metaD$experiment_info)) empty=TRUE
  grpInp = groupsInput()
  if(is_empty(grpInp)) empty=TRUE
  if(any(lapply(grpInp,is_empty) %>% unlist())) empty=TRUE
  rs <- rsUsed()
  req(rs)
  relativeToCol = input$relativeTo
  numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
  expInfo = metaD$experiment_info[,numCols, drop=F]
  rownames(expInfo) = metaD$experiment_info$sample_id
  outputList = rv$outputList

  if(length(rs)<1) empty = TRUE
  if(empty){
    output$ttestResult = renderPrint("no data")
    return(NULL)
  }
  # browser()
  if(relativeToCol == "none"){
    rSums = rep(1,nrow(clusterPatientTable))
  } else {
    if (relativeToCol %in% colnames(metaD$experiment_info)){
      rSums = expInfo[rownames(clusterPatientTable),relativeToCol]
    } else {
      if(relativeToCol %in% names(outputList)){
        rSums = rowSums(clusterPatientTable[,outputList[[relativeToCol]],drop=F])
      } else{
        cat(file = stderr(), "ERROR\n")
        return(NULL)
      }
    }
  }
  names(rSums) = rownames(clusterPatientTable)
  cD = colData(sce)
  cD=cD[cD$cluster_id %in% rs,]
  cellCounts = cD %>% as_tibble() %>% group_by(sample_id) %>% count()
  # save(file = "/pasteur/appa/scratch/bernd/dev.RData", list = c("empty", "grpInp", "rs", "cD", "cellCounts", "relativeToCol", "rSums"))
  # cp =load(file = "/pasteur/appa/scratch/bernd/dev.RData")
  x = cellCounts %>% filter(sample_id %in% grpInp$group1) %>% ungroup()
  x$rsums = rSums[x$sample_id]
  x$val = x$n / x$rsums
  x = x %>% pull(val)
  y = cellCounts %>% filter(sample_id %in% grpInp$group2) %>% ungroup()
  y$rsums = rSums[y$sample_id]
  y$val = y$n / y$rsums
  y = y %>% pull(val)
  req(x); req(y)
  if(sd(x)==0 & sd(y)==0) {
    output$ttestResult = renderPrint("not enough data")
    return(NULL)
  }
  tt = t.test(x,y)
  # if(is.na(tt$p.value)) next()
  # if(tt$p.value < minNodes[[rIdx]][["p.val"]]){
  #   minNodes[[rIdx]][["p.val"]] = tt$p.value
  #   minNodes[[rIdx]][["diff"]] = mean(x) - mean(y)
  #   minNodes[[rIdx]][["node"]] = node
  # }
  output$ttestResult = renderPrint(tt)
})

# updatedoutputList function ----
updatedoutputList <- function(){
  ol <- rv$outputList
  updateSelectInput(session = session, "clusterNameRM", choices = names(ol))
  oldVal = isolate(input$clusterNameSelect)
  updateSelectInput(session = session, "clusterNameSelect", choices = names(ol), selected = oldVal)
  oldval = input$compareStatsTo
  # expInfo = metaD$experiment_info
  numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
  expInfo = metaD$experiment_info[,numCols, drop=F]
  # expInfo = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]
  eI = apply(expInfo,2,as.numeric)
  # browser()
  if(is.na(eI) %>% any()) {
    cat(file = stderr(), "some NAs produced 1")
    browser()
  }
  choices = c("none",colnames(eI), names(ol))
  updateSelectInput(session = session, "compareStatsTo", choices = choices, selected = oldval)
  oldval = input$relativeTo
  choices = c("none",colnames(eI), names(ol))
  updateSelectInput(session = session, "relativeTo", choices = choices, selected = oldval)
  oldval = input$upsetSelection
  updateSelectInput(session = session, "upsetSelection", choices = names(ol), selected = oldval)
  oldval = input$colorbyGroups
  updateSelectInput(session = session, "colorbyGroups", choices = names(ol), selected = oldval)
  oldval = input$groupRM
  updateSelectInput(session = session, "groupRM", choices = names(ol), selected = oldval)

}

# observe applyclusterNumbers ---
shiny::observeEvent(input$applyclusterNumbers, {
  outputList = rv$outputList
  # isolate(input$clusterNumbers)
  updatedoutputList()
})

# obs rmGrp ----
shiny::observeEvent(input$rmGrp,{
  cat(file = stderr(), "obs rmGrp ----\n")
  outputList = rv$outputList
  cName = input$clusterNameRM
  cList = inputClusterNumber()
  if(cName %in% names(outputList)){
    outputList[[cName]] = NULL
    outputList[["Rest"]] = c()
    outputList[["Rest"]] = as.integer(levels(sce$cluster_id)[!levels(sce$cluster_id) %in% unique(unlist(outputList))])
    for(na in names(outputList)){
      if(length(outputList[[na]])==0) outputList[[na]] = NULL
    }
    # outputList = append(outputList, list(cList ))
    # names(outputList)[length(outputList)] = cName
    rv$outputList = outputList
    updatedoutputList()
  }
})


# obs applyName ----
shiny::observeEvent(input$applyName,{
  cat(file = stderr(), "obs applyName ----\n")
  outputList = rv$outputList
  cName = input$clusterName
  cList = inputClusterNumber()
  outputList = rv$outputList
  outputList[[cName]] = cList
  outputList[["Rest"]] = c()
  outputList[["Rest"]] = as.integer(levels(sce$cluster_id)[!levels(sce$cluster_id) %in% unique(unlist(outputList))])
  for(na in names(outputList)){
    if(length(outputList[[na]])==0) outputList[[na]] = NULL
  }
  # browser()
  # outputList = append(outputList, list(cList ))
  # names(outputList)[length(outputList)] = cName
  rv$outputList = outputList
  updatedoutputList()
  currentSelection = input$colorbyGroups
  updateSelectInput(inputId = "colorbyGroups",
                    choices = names(outputList),
                    selected = c(currentSelection,cName)
  )
})

shiny::observeEvent(input$rmGroups,{
  cat(file = stderr(), "obs rmGroups ----\n")
  outputList = rv$outputList
  rmGroups = input$groupRM
  rmCluster = outputList[rmGroups] %>% unlist() %>% unique()
  rs = rsUsed_d() %>% isolate()
  rsUsed(setdiff(rs, rmCluster))
})

# output dend ----
output$dend <- renderPlot({
  dendPlot() %>%
    plot(main = "dendrogram")
})

dendPlot <- reactive({
  cat(file = stderr(), "dend\n")
  rs = rsUsed()
  req(rs)
  selectedPoints()
  labCol = rep("blue",length(labels(dend)))
  labCol[which(labels(dend) %in% (rs))] = "red"
  dend %>% dendextend::set("leaves_col", labCol) %>%
    dendextend::set("leaves_pch", 15) %>%
    dendextend::set("leaves_cex", 1) %>%
    dendextend::set("labels_cex", 0.5)
})

# observe({
#   d <- event_data("plotly_selected")
#   # d2 <- event_data("plotly_selected")
#   if(is.null(d)){return(NULL)}
#   if(nrow(d)==0){return(NULL)}
#   d = d[d$curveNumber==0,]
#   rsUsed(d$pointNumber+1)
# })

# output str ----
output$str <- renderPrint({
  str(dendTable[dendTable$child==input$node[[length(input$node)]],])
})

# observeEvent(input$d2Axes1 ---
shiny::observeEvent(input$d2Axes1,{
  # browser()
  activePlot(as.numeric(input$d2Axes1))
})

shiny::observeEvent(activePlot(),{
  cat(file = stderr(), "active plot: ", activePlot(), "\n")
})

# inputSelect = function ----
inputSelect = function(d, rs, mode){
  req(rs)
  req(d)
  req(d$curveNumber)
  # browser()
  # d = d[d$curveNumber==0,]
  d = d$customdata
  # "view", , "add"

  d = switch(EXPR = mode,
             "remove others" = intersect(d, rs),
             "add" = unique(c(d, rs)),
             "remove" = rs[!rs %in% d],
             d
  )
  d
}
# inputClusterNumber <- reactive ----
inputClusterNumber <- reactive({
  cat(file = stderr(), "inputClusterNumber: ", paste(input$clusterNumbers, collapse = ", "),"\n")
  str_split(input$clusterNumbers, ",")[[1]] %>% as.integer()
}) %>% debounce(1000)

# violinPlotSelection <- reactive ----
violinPlotSelection <- reactive({
  input$violinSelection
})%>% debounce(1000)

# inputClusterNumber() ----
shiny::observe({
  ic = inputClusterNumber()
  # in case some clusters are missing
  ic = as.integer(intersect(ic, colnames(clusterPatientTable)))
  isolate(rsUsed( ic))
})

# rsUsed_d
rsUsed_d <- rsUsed %>% debounce(10000)

# updateTextInput clusterNumbers ----
shiny::observe({
  rs= rsUsed() %>% sort
  outputList = rv$outputList
  updateOL = FALSE
  if(!"selected" %in% names(outputList)){
    updateOL = TRUE
  }
  outputList$selected = rs
  rv$outputList = outputList
  if(updateOL){
    updatedoutputList()
  }
  updateTextInput(inputId = "clusterNumbers", value =  paste(rs, collapse = ", ") )
})


shiny::observeEvent(event_data("plotly_selected", source = "somGrid"),{
  message("somGrid touched")
  # browser()
  rs = rsUsed_d() %>% isolate()
  req(rs)
  d <- event_data("plotly_selected", source = "somGrid")
  if(is.null(d)){return(NULL)}
  d = inputSelect(d, rs, isolate(input$selectMode))
  isolate(rsUsed(d))
})


lapply(seq_len(nPlots), function(i) {
  observeEvent(event_data("plotly_selected", source = paste0("somData", i)), {
    message("som", i, " touched")
    activePlot(i)
    rs <- isolate(rsUsed_d())
    req(rs)
    d <- event_data("plotly_selected", source = paste0("somData", i))
    if (is.null(d)) return(NULL)
    d <- inputSelect(d, rs, isolate(input$selectMode))
    isolate(rsUsed(d))
  })
})
#
# # observeEvent(event_data("plotly_selected", source = "somData1"),{ ----
# #2DO: need to be adjusted for variable list of plots
# shiny::observeEvent(event_data("plotly_selected", source = "somData1"),{
#   message("som1 touched")
#   # browser()
#   activePlot(1)
#   rs = rsUsed_d() %>% isolate()
#   req(rs)
#   d <- event_data("plotly_selected", source = "somData1")
#   if(is.null(d)){return(NULL)}
#   d = inputSelect(d, rs, isolate(input$selectMode))
#   isolate(rsUsed(d))
# })
# shiny::observeEvent(event_data("plotly_selected", source = "somData2"),{
#   message("som2 touched")
#   activePlot(2)
#   rs = rsUsed_d() %>% isolate()
#   req(rs)
#   d <- event_data("plotly_selected", source = "somData2")
#   if(is.null(d)){return(NULL)}
#   d = inputSelect(d, rs, isolate(input$selectMode))
#   # selectedPoints(d)
#   isolate(rsUsed(d))
# })
# shiny::observeEvent(event_data("plotly_selected", source = "somData3"),{
#   message("som3 touched")
#   activePlot(3)
#   rs = rsUsed_d() %>% isolate()
#   req(rs)
#   d <- event_data("plotly_selected", source = "somData3")
#   if(is.null(d)){return(NULL)}
#   d = inputSelect(d, rs, isolate(input$selectMode))
#   # selectedPoints(d)
#   isolate(rsUsed(d))
# })
# shiny::observeEvent(event_data("plotly_selected", source = "somData4"),{
#   message("som4 touched")
#   activePlot(4)
#   rs = rsUsed_d() %>% isolate()
#   req(rs)
#   d <- event_data("plotly_selected", source = "somData4")
#   if(is.null(d)){return(NULL)}
#   d = inputSelect(d, rs, isolate(input$selectMode))
#   # selectedPoints(d)
#   isolate(rsUsed(d))
# })
# shiny::observeEvent(event_data("plotly_selected", source = "somData5"),{
#   message("som5 touched")
#   activePlot(5)
#   rs = rsUsed_d() %>% isolate()
#   req(rs)
#   d <- event_data("plotly_selected", source = "somData5")
#   if(is.null(d)){return(NULL)}
#   d = inputSelect(d, rs, isolate(input$selectMode))
#   isolate(rsUsed(d))
# })
# shiny::observeEvent(event_data("plotly_selected", source = "somData6"),{
#   # activePlot(6)
#   message("som6 touched")
#   cat(file = stderr(), "obs: somData6\n")
#   rs = rsUsed_d() %>% isolate()
#   req(rs)
#   d <- event_data("plotly_selected", source = "somData6")
#   if(is.null(d)){return(NULL)}
#   d = inputSelect(d, rs, isolate(input$selectMode))
#   isolate(rsUsed(d))
# })
shiny::observeEvent(event_data("plotly_selected", source = "tsne"),{
  message("tsne touched")
  rs = rsUsed_d() %>% isolate()
  req(rs)
  d <- event_data("plotly_selected", source = "tsne")
  if(is.null(d)){return(NULL)}
  # browser()
  # save(file = "/pasteur/appa/scratch/bernd/event.RData", list = c("d", "rs"))
  #
  d = inputSelect(d, rs, isolate(input$selectMode))
  isolate(rsUsed(d))
})
shiny::observeEvent(event_data("plotly_selected", source = "umap"),{
  message("umap touched")
  rs = rsUsed_d() %>% isolate()
  req(rs)
  d <- event_data("plotly_selected", source = "umap")
  if(is.null(d)){return(NULL)}
  d = inputSelect(d, rs, isolate(input$selectMode))
  isolate(rsUsed(d))
})
shiny::observeEvent(event_data("plotly_selected", source = "pca"),{
  message("pca touched")
  rs = rsUsed_d() %>% isolate()
  req(rs)
  d <- event_data("plotly_selected", source = "pca")
  if(is.null(d)){return(NULL)}
  d = inputSelect(d, rs, isolate(input$selectMode))
  isolate(rsUsed(d))
})

# observe input$node dendTable ----
shiny::observe({
  cat(file = stderr(), "obs rs\n")
  rs = labels(eval(parse(text=paste0("dend",dendTable[dendTable$child==input$node[[length(input$node)]],"indexString"]))))
  rsUsed(rs) %>% isolate()
})


# output$selected  ----
output$selected <- renderPrint({
  cat(file = stderr(), "selected\n")

  rs <- rsUsed()
  req(rs)
  rs
})

# cellCounts ----
output$cellCounts <- DT::renderDT(options = list(lengthChange = FALSE,
                                                 scrollX = TRUE),
                                  {
                                    cat(file = stderr(), "cellCounts\n")
                                    cst = input$compareStatsTo
                                    sN = input$singleNode %>% as.integer()
                                    rs <- rsUsed()
                                    req(rs)

                                    outputList = rv$outputList
                                    print("cellCounts\n")

                                    # browser()
                                    req(clusterPatientTable)
                                    rSums = rowSums(clusterPatientTable[,rs,drop=F])
                                    names(rSums) = rownames(clusterPatientTable)
                                    rSums = as.data.table(t(rSums))
                                    if(!cst=="none"){
                                      if(cst %in% colnames(metaD$experiment_info)){
                                        # expInfo = metaD$experiment_info
                                        numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
                                        expInfo = metaD$experiment_info[,numCols, drop=F]
                                        eI = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]
                                        eI = apply(eI,2,as.numeric) %>% as.data.frame()
                                        rownames(eI) = metaD$experiment_info$sample_id
                                        if(is.na(eI) %>% any()) {
                                          cat(file = stderr(), "some NAs produced 2\n")
                                          browser()
                                        }
                                        # ctStats = eI[,cst]
                                        ctStats =  as.data.table(t(eI[colnames(rSums),cst]))
                                        colnames(ctStats) = colnames(rSums)
                                        rSums = data.table::rbindlist(list(rSums, ctStats))
                                        # rbind(rSums, eI[colnames(rSums),cst])
                                        rownames(rSums) = c("selection", cst)
                                      } else if(cst %in% names(outputList)){
                                        ctStats =  as.data.table(t(rowSums(clusterPatientTable[,outputList[[cst]],drop=F])))
                                        colnames(ctStats) = rownames(clusterPatientTable)
                                        rSums = data.table::rbindlist(list(rSums, ctStats))
                                        # rSums = data.table::rbindlist(list(rSums, as.data.table(t())))
                                        rownames(rSums) = c("selection", cst)
                                      } else{
                                        cat(file = stderr(), "should not happen\n")
                                        # browser()
                                      }
                                    }
                                    rn = rownames(rSums)
                                    ctStats = as.data.table(t(rowSums(clusterPatientTable[,sN,drop=F])))
                                    colnames(ctStats) = rownames(clusterPatientTable)
                                    rSums = data.table::rbindlist(list(rSums, ctStats))
                                    # rbind(rSums, rowSums(clusterPatientTable[,sN,drop=F]))
                                    rownames(rSums) = c(rn, paste("SOM node", sN))
                                    # save(file = "/pasteur/appa/scratch/bernd/Rtest3.Rdata",
                                    #      list = c("clusterPatientTable", "rs", "cst", "rSums", "sN"))
                                    # cp =load(file = "/pasteur/appa/scratch/bernd/Rtest3.Rdata")
                                    print(length(rSums))
                                    print(class(rSums))
                                    print(nrow(rSums))
                                    # cat(rSums)
                                    rSums
                                  })

# countBar ----
output$CountBar <-renderPlot({
  countBarPlot()
})

countBarPlot <- reactive({
  cat(file = stderr(), "countBar\n")
  cst = input$compareStatsTo
  # browser()
  outputList = rv$outputList
  rs <- rsUsed()
  req(rs)
  groupsInput = groupsInput()
  req(clusterPatientTable)
  countBarPlotFunc(rs, clusterPatientTable, cst, sce, outputList, groupsInput)

})

# cellPercentages ----
output$cellPercentages <- renderPrint({
  cat(file = stderr(), "cellPercentages\n")

  outputList = rv$outputList
  rs <- rsUsed()
  req(rs)

  req(clusterPatientTable)
  relativeToCol = input$relativeTo
  expInfo = metadata(sce)$experiment_info
  rownames(expInfo) = expInfo$sample_id
  expInfo = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=F]
  eI = apply(expInfo,2,as.numeric)

  rSums <- compute_relative_counts(
    clusterPatientTable, rs, relativeToCol, expInfo, outputList
  )
  names(rSums) <- rownames(clusterPatientTable)

  # if(relativeToCol == "none"){
  #   rSums = rowSums(clusterPatientTable[,rs,drop=F])/rowSums(clusterPatientTable[,,drop=F])*100
  # } else {
  #   if (relativeToCol %in% colnames(metadata(sce)$experiment_info)){
  #     rSums =rowSums(clusterPatientTable[,rs,drop=F])/expInfo[rownames(clusterPatientTable),relativeToCol]*100
  #   } else {
  #     if(relativeToCol %in% names(outputList)){
  #       rSums = rowSums(clusterPatientTable[,rs,drop=F])/rowSums(clusterPatientTable[,outputList[[relativeToCol]],drop=F])*100
  #     } else{
  #       cat(file = stderr(), "ERROR\n")
  #     }
  #   }
  # }
  # names(rSums) = rownames(clusterPatientTable)
  # rSums$id = rownames(clusterPatientTable)
  # names(rSums) = rownames(clusterPatientTable)
  # print(rSums)
  print(noquote(formatC(signif(rSums, digits=2), digits=2, format="fg", flag="#")))
})

#PercentBar ----
output$PercentBar <-renderPlot({
  PercentBarPlot()
})

PercentBarPlot <- reactive({
  cat(file = stderr(), "PercentBar\n")
  outputList = rv$outputList
  rs <- rsUsed()
  req(rs)
  req(clusterPatientTable)
  relativeToCol = input$relativeTo
  groupsInput = groupsInput()
  PercentBarPlotFunc(sce, relativeToCol, clusterPatientTable, rs, outputList, group, groupsInput)
})


### zooming ----
zoomFunc <- function(zoom,  plotIdx) {
  dimSelectionInternal =  dimSelection()
  req(dimSelectionInternal)
  cat(file = stderr(), "zoomFunc: ", plotIdx, "\n")
  if(all(!is.null(zoom), "xaxis.range[0]" %in% names(zoom) , na.rm = T)) {
    # browser()
    rezoom = FALSE
    if (all( c(dimSelectionInternal[[plotIdx]]$xzoom[1] > zoom$`xaxis.range[0]`,
               !is.null(dimSelectionInternal[[plotIdx]]$xzoom[1])),na.rm = T)) rezoom = TRUE
    if(rezoom){
      # cat(file = stderr(), "rezoom\n")

      # browser()
      dimSelectionInternal[[plotIdx]]$xzoom <- c(NULL,NULL)
      dimSelectionInternal[[plotIdx]]$yzoom <- c(NULL,NULL)
    } else {
      # cat(file = stderr(), "zoom\n")
      # cat(file = stderr(), paste(dimSelection[[plotIdx]]$xzoom[1], zoom$`xaxis.range[0]`))
      # cat(file = stderr(), "\n")
      # browser()
      dimSelectionInternal[[plotIdx]]$xzoom <- c(zoom$`xaxis.range[0]`,zoom$`xaxis.range[1]`)
      dimSelectionInternal[[plotIdx]]$yzoom <- c(zoom$`yaxis.range[0]`,zoom$`yaxis.range[1]`)
    }
    dimSelection(dimSelectionInternal)
    isolate(triggerRedraw(triggerRedraw() + 1))
  }
}

# obs somData ----
shiny::observe({
  cat(file = stderr(), "obs Event somData 1\n")

  plotIdx = 1
  zoom <- event_data("plotly_relayout", source = paste0("somData",plotIdx))
  zoomFunc(zoom, plotIdx)
})
shiny::observe({
  cat(file = stderr(), "obs Event somData 2\n")

  plotIdx = 2
  zoom <- event_data("plotly_relayout", source = paste0("somData",plotIdx))
  zoomFunc(zoom, plotIdx)
})
shiny::observe({
  cat(file = stderr(), "obs Event somData 3\n")

  plotIdx = 3
  zoom <- event_data("plotly_relayout", source = paste0("somData",plotIdx))
  zoomFunc(zoom, plotIdx)
})
shiny::observe({
  cat(file = stderr(), "obs Event somData 4\n")

  plotIdx = 4
  zoom <- event_data("plotly_relayout", source = paste0("somData",plotIdx))
  zoomFunc(zoom, plotIdx)
})
shiny::observe({
  cat(file = stderr(), "obs Event somData 5\n")

  plotIdx = 5
  zoom <- event_data("plotly_relayout", source = paste0("somData",plotIdx))
  zoomFunc(zoom, plotIdx)
})
shiny::observe({
  cat(file = stderr(), "obs Event somData 6\n")

  plotIdx = 6
  zoom <- event_data("plotly_relayout", source = paste0("somData",plotIdx))
  zoomFunc(zoom, plotIdx)
})
# observe samples2plot. ---
sample2PlotDb <- reactive(input$samples2plot) %>% debounce(500)

shiny::observe({
  sampleIds = sample2PlotDb()
  dfPlot <<- data.frame(t(assays(sce_subsampled[,sce_subsampled$sample_id %in% sampleIds])[[1]]))
})

shiny::observeEvent(dimSelection(),{
  dimSelection = dimSelection()
  for(idx in seq_along(dimSelection)){
    cat(file = stderr(), paste(idx, dimSelection[[idx]], collapse = " "))
  }
})

shiny::observeEvent(event_data("plotly_selected", source = "scatterPlot"),{
  cat(file = stderr(), "obs Event scatter\n")

  sc=sce
  rs = isolate(rsUsed_d())
  dimSelection =  dimSelection()
  sampleIds = isolate(input$samples2plot)
  # browser()
  # req(rs)
  d <- event_data("plotly_selected", source = "scatterPlot")
  plotIdx = activePlot()
  if(is.null(d)){return(NULL)}
  req(d$curveNumber)
  d = d[d$curveNumber==1,]
  # d = d$pointNumber+1
  # "view", , "add"
  minx = min(d$x)
  maxx = max(d$x)
  miny = min(d$y)
  maxy = max(d$y)
  ids = which( dfPlot[,dimSelection[[plotIdx]]$dims[1] %>% make.names()] > minx &
                 dfPlot[,dimSelection[[plotIdx]]$dims[1] %>% make.names()]  < maxx &
                 dfPlot[,dimSelection[[plotIdx]]$dims[2] %>% make.names()] > miny &
                 dfPlot[,dimSelection[[plotIdx]]$dims[2] %>% make.names()]  < maxy  )
  ids = colData(sce_subsampled[,sce_subsampled$sample_id %in% sampleIds])[ids,"cluster_id"]
  # activePlot(1)
  # rs = rsUsed()
  # req(rs)
  # d <- event_data("plotly_selected", source = "somData1")
  # if(is.null(d)){return(NULL)}
  ids = switch(EXPR = isolate(input$selectMode),
               "remove others" = intersect(ids, rs),
               "add" = unique(c(ids, rs)),
               "remove" = rs[!rs %in% ids],
               ids
  )
  isolate(rsUsed(ids))
  # isolate(rsUsed(colData(sce)[ids,"cluster_id"] %>% unique()))

  # rsUsed( unique(colData(sce[,colData(sce)$cluster_id %in% rs])$"cluster_id" ))
  # # d
  # browser()
  # # d = inputSelect(d, rs, isolate(input$selectMode))
  # rsUsed(d)
})


###   output$somData* ----
# Shared base plot for each dimSelection slot, cached per channel pair.
# This avoids rebuilding plotSOMScatter from scratch for each of the 6 outputs.
somBasePlots <- lapply(seq_len(nPlots), function(plotIdx) {
  reactive({
    dimSel <- dimSelection()
    req(dimSel)
    dims <- dimSel[[plotIdx]]$dims
    plotSOMScatter(x = sce,
                   chs = c(dims[1], dims[2]),
                   pointSize = "max",
                   color_by = "n",
                   xRN = sceRN, xCN = sceCN) +
      scale_colour_gradientn(colours = viridis::viridis(9))
  }) %>%
    bindCache(dimSelection()[[plotIdx]]$dims)
})

# Cached projection data frames for the showGroups SOM view.
# This avoids repeating ggplot_build() + left_join() every time rs or
# colorbyGroups change; only the final ggplot construction remains.
projectionDfs <- lapply(seq_len(nPlots), function(plotIdx) {
  reactive({
    pp1 <- somBasePlots[[plotIdx]]()
    req(pp1)
    buildProjectionDf(pp1, plotIdx, dimSelection(), sce)
  }) %>%
    bindCache(dimSelection()[[plotIdx]]$dims)
})

lapply(seq_len(nPlots), function(i) {
  local({
    plotIdxLocal <- i
    output[[paste0("somData", plotIdxLocal)]] <- renderPlotly({
      cat(file = stderr(), "som", plotIdxLocal, "\n")
      colorbyGroups <- input$colorbyGroups
      selectedUpdate2()
      showGroups <- input$showGroups
      dimSelection <- dimSelection()
      rs <- rsUsed_d()
      req(rs)
      triggerRedraw()
      plotIdx <- if (plotIdxLocal == 6) activePlot() else plotIdxLocal
      pp1 <- somBasePlots[[plotIdx]]()
      projectionDf <- if (showGroups) projectionDfs[[plotIdx]]() else NULL
      p3 <- somPlot(pp1, plotIdx, rs, colorbyGroups, showGroups,
                    dimSelection = dimSelection, sce = sce, metaD = metaD,
                    outputList = rv$outputList, projectionDf = projectionDf)
      ggplotly(p3, source = paste0("somData", plotIdxLocal), tooltip = "") %>%
        layout(showlegend = F, dragmode = "select") %>%
        config(renderer = "webgl") %>%
        event_register("plotly_selected") %>%
        event_register("plotly_relayout")
    })
  })
})

# color selections ----
# observeEvent(input,{
#   currentSelection = input$colorbyGroups
#   outputList = rv$outputList
#   updateSelectInput(inputId = "colorbyGroups",
#                     choices = names(outputList),
#                     selected = currentSelection
#                     )
# })

dimRedSelection <- reactive({
  retVal = input$dimRedSelection
  retVal
}) %>% debounce(1000)

# tsne ----
tsne <- reactive({
  dimRedCols = dimRedSelection()
  perplexity = input$perplexity
  tsne = tsneFunc(dimRedSelection = dimRedCols, perplexity = perplexity, sce, somCodesName)
  return(tsne)
}) %>%
  debounce(1000) %>%
  bindCache(dimRedSelection(), input$perplexity)

umap <- reactive({
  seed = 1
  dimRedCols = input$dimRedSelection
  set.seed(seed)
  pumap = umap::umap.defaults
  pumap$n_neighbors = input$n_neighbors
  um = umap::umap(metaD[[somCodesName]][,dimRedCols],config = pumap)
  return(um)
}) %>%
  debounce(1000) %>%
  bindCache(input$dimRedSelection, input$n_neighbors)

pca <- reactive({
  seed = 1
  dimRedCols = input$dimRedSelection
  set.seed(seed)
  pca = prcomp(t(metaD[[somCodesName]][,dimRedCols]), scale =F,rank. = 2)
  return(pca)
}) %>%
  debounce(1000) %>%
  bindCache(input$dimRedSelection)

output$tsne <- renderPlotly({
  p3 = tsnePlot()
  showLegend = input$showlegend
  # browser()
  retVal = ggplotly(p3, source = paste0("tsne"), tooltip="text")

  if(showLegend){
    retVal = retVal %>% plotly::layout(legend=list(x=0, y=-3 ,
                                                   xanchor='left',
                                                   yanchor='bottom',
                                                   orientation='h'))
  }else {
    retVal = retVal %>% layout(
      showlegend = F
    )
  }

  retVal = retVal %>%
    layout(dragmode = "select") %>%
    event_register("plotly_selected") %>%
    event_register("plotly_relayout")
  retVal
})

tsnePlot <- reactive({
  cat(file = stderr(), "tsne\n")
  selectedUpdate2()
  rs = rsUsed_d()
  req(rs)
  triggerRedraw()
  tsne = tsne()
  # dimSelection =  dimSelection()

  df = as.data.frame(tsne$Y)
  names(df) = c("tsne1", "tsne2")
  p3 = drawProjection(df, rs, colorbyGroups = input$colorbyGroups, sce = sce, outputList = outputList)
  return(p3)
})

output$umap <- renderPlotly({
  p3 = umapPlot()
  showLegend = input$showlegend
  retVal = ggplotly(p3, source = paste0("umap"), tooltip="text")
  if(showLegend){
    retVal = retVal %>% plotly::layout(legend=list(x=0, y=-3 ,
                                                   xanchor='left',
                                                   yanchor='bottom',
                                                   orientation='h'))
  }else {
    retVal = retVal %>% layout(
      showlegend = F
    )
  }

  retVal = retVal %>%
    layout(dragmode = "select") %>%
    event_register("plotly_selected") %>%
    event_register("plotly_relayout")
  retVal
})


umapPlot <- reactive({
  cat(file = stderr(), "umap\n")
  umap <- umap()
  # plot(um$layout)
  selectedUpdate2()
  rs = rsUsed_d()
  req(rs)
  triggerRedraw()
  df = as.data.frame(umap$layout)
  names(df) = c("umap1", "umap2")

  p3 = drawProjection(df, rs, colorbyGroups = input$colorbyGroups, sce = sce, outputList = outputList)
  return(p3)
})

output$pca <- renderPlotly({
  pca <- pca()
  rs = rsUsed_d()
  req(rs)
  triggerRedraw()
  req(pca)
  df = as.data.frame(pca$rotation)
  colnames(df) = c("pc1", "pc2")
  showLegend = input$showlegend

  p3 = pcaPlot()
  retVal = ggplotly(p3, source = paste0("pca"), tooltip="text")

  if(showLegend){
    retVal = retVal %>% plotly::layout(legend=list(x=0, y=-3 ,
                                                   xanchor='left',
                                                   yanchor='bottom',
                                                   orientation='h'))
  }else {
    retVal = retVal %>% layout(
      showlegend = F
    )
  }
  retVal = retVal %>%
    layout(dragmode = "select") %>%
    event_register("plotly_selected") %>%
    event_register("plotly_relayout")
  return(retVal)
  # retVal %>%
  #   add_trace(x=df$pc1, y=df$pc2, text=~paste("cluster: ", df$cluster),
  #             hoverinfo = 'text', mode = "none")
})

pcaPlot <- reactive({
  cat(file = stderr(), "pca\n")
  pca <- pca()
  selectedUpdate2()
  rs = rsUsed_d()
  req(rs)
  triggerRedraw()
  req(pca)
  df = as.data.frame(pca$rotation)
  colnames(df) = c("pc1", "pc2")
  p3 = drawProjection(df, rs, colorbyGroups = input$colorbyGroups, sce = sce, outputList = outputList)

  return(p3)
})


output$somClusters <- renderText({
  cat(file = stderr(), "somClusters\n")

  rs = rsUsed()
  req(rs)
  rs = as.integer(intersect(rs, colnames(clusterPatientTable)))
  # the row it should be plotted
  as.integer(rs) / 40+1
  rowElement=list()
  lapply(40:1, function(x){rowElement[[x]] = paste(rs[as.integer(rs / 40)+1 == x],collapse = ", ")}) %>% unlist() %>% paste(collapse = "\n")
  # paste(sort(rs,decreasing = T), collapse = ", ")
})


### output$scatter ----
output$scatter <- renderPlotly({
  rs = rsUsed()
  req(rs)
  dimSelection =  dimSelection()
  sampleIds = input$samples2plot
  # browser()
  plotIdx = activePlot()
  cidIdx = colData(sce_subsampled)$cluster_id %in% rs & colData(sce_subsampled)$sample_id %in% sampleIds
  if(length(cidIdx)<1)return(NULL)

  pp = scatterPlot()
  # here we add a grid for selecting points
  minx = df[cidIdx,dimSelection[[plotIdx]]$dims[1] %>% make.names()] %>% min()
  maxx =  df[cidIdx,dimSelection[[plotIdx]]$dims[1] %>% make.names()] %>% max()
  xpoints = seq(from = minx, to=maxx, length=100) %>% rep(100) %>% sort()
  miny = df[cidIdx,dimSelection[[plotIdx]]$dims[2] %>% make.names()] %>% min()
  maxy =  df[cidIdx,dimSelection[[plotIdx]]$dims[2] %>% make.names()] %>% max()
  ypoints = seq(from = miny, to=maxy, length=100)%>% rep(100)
  # browser()
  pp = ggplotly(pp, source = "scatterPlot") %>% add_trace(
    x = xpoints,
    y = ypoints,
    mode='markers',
    type = 'scatter',
    fill = 'none',
    fillcolor = '#e763fa',
    opacity = 0.01 #,     # size=1.9,

  ) %>%
    layout(dragmode = "select") %>%
    event_register("plotly_selected") %>%
    event_register("plotly_relayout")
  return(pp)
})

scatterPlot <- reactive({
  cat(file = stderr(), "scatter\n")

  rs = rsUsed()
  req(rs)
  dimSelection =  dimSelection()
  sampleIds = input$samples2plot
  browser()
  plotIdx = activePlot()
  cidIdx = colData(sce_subsampled)$cluster_id %in% rs & colData(sce_subsampled)$sample_id %in% sampleIds
  if(length(cidIdx)<1)return(NULL)
  pp = plotScatterBJ(rowNames = sce_subsampledRN,
                     x = sce_subsampled[,cidIdx],
                     chs = c(dimSelection[[plotIdx]]$dims[1], dimSelection[[plotIdx]]$dims[2]),
                     bins = 200) +
    xlim(c((dimSelection[[plotIdx]]$xlim[1] %>% as.numeric() %>% sign()) * (dimSelection[[plotIdx]]$xlim[1] %>% as.numeric() %>% abs() %>% ceiling()),
           (dimSelection[[plotIdx]]$xlim[2] %>% as.numeric() %>% sign()) * (dimSelection[[plotIdx]]$xlim[2] %>% as.numeric() %>% ceiling()))) +
    ylim(c((dimSelection[[plotIdx]]$ylim[1] %>% as.numeric() %>% sign()) * (dimSelection[[plotIdx]]$ylim[1] %>% as.numeric() %>% abs() %>% ceiling()),
           (dimSelection[[plotIdx]]$ylim[2] %>% as.numeric() %>% sign()) * (dimSelection[[plotIdx]]$ylim[2] %>% as.numeric() %>% ceiling())))

  return(pp)
})


# somRaster ----
output$somRaster = renderPlot({
  message("somRasterPlot(): ", class(somRasterPlot()))
  res = somRasterPlot()
  raster::plot(res[[1]], addfun = res[[2]], maxnl=80)
})

# somRasterSelect ----
output$somRasterSelect = renderPlotly({
  cat(file = stderr(), "som raster\n")

  rs = rsUsed()
  req(rs)

  data.points = expand.grid(seq(somRasterObj@nrows), seq(somRasterObj@ncols))
  colnames(data.points) = c("x", "y")
  p3 = ggplot(data.points, aes(x,y, customdata=1:nrow(data.points))) + geom_point() +
    geom_point(data=data.points[rs,],
               aes(x=x,
                   y=y, customdata=rs),
               color='red',
               size=0.9)
  ggplotly(p3, source = "somGrid") %>%
    layout(showlegend = F) %>%
    layout(dragmode = "select") %>%
    event_register("plotly_selected") %>%
    event_register("plotly_relayout")

})

somRasterPlot = reactive({
  cat(file = stderr(), "som raster\n")

  rs = rsUsed()
  req(rs)
  xy = somRasterData[rs,c("x","y")]
  fun <- function() {
    points(xy, cex=2)
    points(xy, pch=3, col='red')
  }
  # p = ggplotify::as.ggplot(function() plot(somRasterObj, addfun = fun))
  # p
  return(list(somRasterObj, addfun = fun))
})

## VlnPlot ----

output$VlnPlot = renderPlot({
  vlnPlot()
})
vlnPlot <- reactive({
  cat(file = stderr(), "VlnPlot\n")
  req(input$violinBox)  # only compute while violin box is expanded

  # observe
  # this changes the groups
  input$applyName
  input$rmGrp
  violinSelection = violinPlotSelection()
  upsetSelection = input$upsetSelection
  outputList = rv$outputList
  req(outputList)
  # save(file = "vln.Rdata", list=ls())
  ##                  Cell Idents  Feat     Expr
  # Idents = named groups of clusters
  # Feat = marker
  # Expr = SOM value
  # cell = SOM node
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
    wide$somNode = factor(outputList[[na]], levels = 1:nrow(sce@metadata[[somCodesName]]))
    long = gather(wide, marker,expr,-somNode, factor_key=TRUE)
    long$grpName = na
    data = rbind(data, long)
  }
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
  p
})
## VlnPlot2 ----
output$VlnPlot2 = renderPlot({
  VlnPlot2()
})


VlnPlot2 = reactive({
  cat(file = stderr(), "VlnPlot2\n")
  req(input$violinBox)  # only compute while violin box is expanded
  # observe
  # this changes the groups
  input$applyName
  input$rmGrp
  upsetSelection = input$upsetSelection
  violinSelection = violinPlotSelection()
  outputList = rv$outputList
  req(outputList)
  # save(file = "vln.Rdata", list=ls())
  ##                  Cell Idents  Feat     Expr
  # Idents = named groups of clusters
  # Feat = marker
  # Expr = SOM value
  # cell = SOM node
  p = plotViolin2Func(sce, somCodesName, violinSelection, upsetSelection, outputList)
  return(p)
})

## close app ----
shiny::observeEvent(input$close, {
  shinyjs::js$closeWindow()
  stopApp()
})

## upset plot ----
output$UpSet = renderPlot({
  upSetPlot()
})

selectedUpdate <- reactiveVal(value=0)

shiny::observe({
  rs = rsUsed_d()
  if("selected" %in% input$upsetSelection){
    isolate(selectedUpdate(selectedUpdate() +1))
  }
})

upSetPlot = reactive({
  cat(file = stderr(), "UpSet: ","\n")
  req(input$upsetBox)  # only compute while UpSet box is expanded
  # inputList = reactiveValuesToList(input)
  # save(file = "UpSet.RData", list = c(ls()))
  # cp = load("UpSet.RData")
  selectedUpdate()
  input$applyName
  upsetSelection = input$upsetSelection
  input$rmGrp
  outputList = rv$outputList
  upsetPlotFunc(upsetSelection, outputList, sce)


})


output$downloadPlots <- downloadHandler(
  filename = function() {
    req(input$clusterNameSelect)
    paste(input$clusterNameSelect, '.pdf', sep='')
  },
  content = function(file) {
    dimSelection =  dimSelection()
    rs = c(1,2,3,4,8)
    rs = rsUsed()
    req(rs)
    cat("class:", class(upSetPlot()))
    pdf(file = file)

    plotIdx = 1

    message("printing dendPlot", class(dendPlot()))
    dendPlot() %>%
      plot(main = "dendrogram")
    message("printing countBarPlot")
    print(countBarPlot())

    message("printing PercentBarPlot")
    print(PercentBarPlot())


    for(plotIdx in seq(nPlots)){
      message("printing somScatter", plotIdx)
      pp1 = plotSOMScatter(x=sce,
                           chs=c(dimSelection[[plotIdx]]$dims[1], dimSelection[[plotIdx]]$dims[2]),
                           pointSize = "max",
                           color_by ="n",xRN = sceRN, xCN = sceCN ) +
        scale_colour_gradientn(colours=viridis::viridis(9))
      print(ggsomPlot(pp1, plotIdx, rs, dimSelection, sce=sce, metaD = metaD))
    }
    message("printing tsnePlot")
    print(tsnePlot())
    message("printing umapPlot")
    print(umapPlot())
    message("printing pcaPlot")
    print(pcaPlot())

    message("printing scatterPlot")
    print(scatterPlot())
    message("printing somRasterPlot")
    res = somRasterPlot()
    raster::plot(res[[1]], addfun = res[[2]], maxnl=80)

    message("printing vlnPlot")
    print(vlnPlot())

    message("printing VlnPlot2")
    print(VlnPlot2())
    message("printing upSetPlot")
    print(upSetPlot())
    dev.off()
    message("done printing")
  }
)

# observer zoom dimSelection ----
dimSelection = reactiveVal(list())
shiny::observe({
  cat(file = stderr(), "\n---changedimSelection\n\n")
  dimSelection = list()
  for(idx in 1:nPlots){
    d1 <- input[[paste0("d",idx,".1")]]
    d2 <- input[[paste0("d",idx,".2")]]
    lim1 <- channelLimits[[d1]]
    lim2 <- channelLimits[[d2]]
    dimSelection[[idx]] = list(
      dims = c(d1, d2),
      xlim = c(lim1["min"], lim1["max"]),
      ylim = c(lim2["min"], lim2["max"]),
      xzoom = c(NULL,NULL),
      yzoom = c(NULL,NULL)
    )
  }
  dimSelection(dimSelection)
})

# output$scatterPoints <- renderPlotly({
#   if(!input$showPoints) return(NULL)
#   rs = rsUsed()
#   req(rs)
#   # browser()
#   plotIdx = activePlot()
#   # browser()
#
#   pp = ggplot( data = df[colData(sce)$cluster_id %in% rs,],
#                mapping = aes_string(x=dimSelection[[plotIdx]]$dims[1] %>% make.names(),
#                                     y=dimSelection[[plotIdx]]$dims[2]%>% make.names())) + geom_point()
#
#   ggplotly(pp, source = "scatterPlot") %>%
#     layout(dragmode = "select") %>%
#     event_register("plotly_selected") %>%
#     event_register("plotly_relayout")
#
# })
    }



  }

  # sce = sce# main input has to contain:
  #                      sce_subsampled = sce_subsampled # subsampled sce object
  #                      outputList = outputList # list of named nodes
  #                      colTree = NULL # Tree object to plot
  #                      dList = dList
  #                      dend = NULL
  #                      dendTable = NULL
  #                      clusterPatientTable = clusterPatientTable
  #                      somCodesName = "SOM_codes" # SOM_codes.1
  #                      nPlots = 6
  #                      somRasterData = somRasterData
  #                      somRasterObj = somRasterObj
  #                      # env = environment()
  #
  #
  #
  #   shinyApp(ui = ui, server = server)

  return(list(ui, server))
}



  