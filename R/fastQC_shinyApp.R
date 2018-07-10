#' Run fastQC shiny app
#'
#' @description Returns a shiny app interface to parse many fastQC objects
#'
#' @details Currently some plots can take a while to render if the \code{FastqcDataList} passed to
#' \code{fastqcInput} has many elements
#'
#' @param fastqcInput can be a \code{FastqcFileList}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files.
#'
#'
#' @return UI data for fastQC shiny.
#'
#' @import shinydashboard
#' @importFrom reshape2 dcast
#' @importFrom magrittr %>%
#' @importFrom plotly layout
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom plotly event_data
#' @import shiny
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles shinyDirChoose
#' @importFrom shinyFiles shinyDirButton
#' @importFrom shinyFiles parseFilePaths
#' @importFrom shinyFiles parseDirPath
#' @importFrom shinyFiles getVolumes
#'
#' @examples
#' \dontrun{
#' #' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # Run the Shiny app
#' fastqcShiny(fdl)
#' }
#'
#' @export
#' @rdname fastqcShiny
fastqcShiny <- function(fastqcInput = NULL){
  
  if(!is.null(fastqcInput)){
    stopifnot(length(fastqcInput) > 1)
    stopifnot(class(fastqcInput) == "character" | class(fastqcInput)[1] == "FastqcDataList")
  }
  
  menuItemLogic <- function(flags){
    
    menuLogic <- list()
    if(all(flags$Status == "PASS")){
      
      menuLogic[[1]] <- "PASS"
      menuLogic[[2]] <- "green"
    }
    else {
      menuLogic[[1]] <- "WARN"
      menuLogic[[2]] <- "yellow"
    }
    
    if(all(flags$Status == "FAIL")){
      menuLogic[[1]] <- "FAIL"
      menuLogic[[2]] <- "red"
    } 

    # Fail values
    menuLogic[[3]] <- c(sum(flags$Status == "FAIL"), length(flags$Status))

    # Warn values
    menuLogic[[4]] <- c(sum(flags$Status == "WARN"), length(flags$Status))

    # Pass values
    menuLogic[[5]] <- c(sum(flags$Status == "PASS"), length(flags$Status))

    
    
    menuLogic
  }
  
  
  renderValBox <- function(count, status, ic, c){
    renderValueBox({
      valueBox(
        value = paste(count[1], count[2], sep = "/"), subtitle = status, icon = icon(ic,
                                                                                     class = "fa-lg"),
        color = c
      )
    })
  }
  
  
  
  body <- dashboardBody(
    tabItems(tabItem(tabName = "BS",
                     column(width = 2,
                            box(h5("Choose FastQC Report:"),
                                shinyFiles::shinyFilesButton(id = "files", label = "Choose files", multiple = TRUE, title = ""),
                                br(),
                                textOutput("report"),
                                br(),
                                checkboxInput("Sumcluster", "Cluster", value = TRUE), 
                                collapsible = TRUE, width = NULL, title = "Options")),
                     box(h1("Summary of fastQC Flags"),
                         h5("Heatmap of fastQC flags (pass, warning or fail) for each fastQC report"),
                         plotlyOutput("SummaryFlags"), width = 10)),
             tabItem(tabName = "TS",
                     box(checkboxInput("showDup", "Show Duplicated?", value = FALSE), 
                         collapsible = TRUE, width = 2, title = "Options"),
                     box(h1("Total Sequences"),
                         h5("Total number of unique and duplicated reads in each sample"),
                         plotlyOutput("ReadTotals"), width = 10)
             ),
             tabItem(tabName = "BQ",
                     column(width = 2, box(radioButtons(inputId="BQplotValue", label="Base Quality",
                                                        choices=c("Mean","Median"), selected = "Mean"),
                                           checkboxInput("BQcluster", "Cluster", value = TRUE),
                                           collapsible = TRUE, width = NULL, title = "Options"),
                     valueBoxOutput("BQboxP",width = NULL),
                     valueBoxOutput("BQboxW",width = NULL),
                     valueBoxOutput("BQboxF",width = NULL)),
                     box(h1("Per Base Sequence Quality"),
                         h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
                         h5("Click sidebar on heatmap to change line plots"),
                         plotlyOutput("baseQualHeatmap"),
                         br(),
                         plotlyOutput("BaseQualitiesSingle"), width = 10)
             ),
             tabItem(tabName = "SQ",
                     column(width = 2, box(
                       radioButtons(inputId="SQType", label="Sequence Quality",
                                    choices=c("Frequency","Counts"), selected = "Frequency"),
                       checkboxInput("SQcluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ), 
                     valueBoxOutput("SQboxP", width = NULL),
                     valueBoxOutput("SQboxW", width = NULL),
                     valueBoxOutput("SQboxF", width = NULL)),
                     box(
                       h1("Per Sequence Quality Scores"),
                       h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
                       h5("Click sidebar on heatmap to change line plots"),
                       plotlyOutput("seqQualHeatmap"),
                       br(),
                       plotlyOutput("SeqQualitiesSingle"), width = 10
                     )
             ),
             tabItem(tabName = "SC",
                     column(width = 2, box(
                       checkboxInput("SCcluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ),
                     valueBoxOutput("SCboxP", width = NULL),
                     valueBoxOutput("SCboxW", width = NULL),
                     valueBoxOutput("SCboxF", width = NULL)),
                     box(
                       h1("Per Base Sequence Content"),
                       h5("Per base sequence content in each sample, colours at each base indicate sequence bias"),
                       h5("G = Black, A = Green, T = Red, C = Blue"),
                       plotlyOutput("SCHeatmap"),
                       br(),
                       plotlyOutput("SCsingle"),
                       width = 10
                     )
             ),
             tabItem(tabName = "GC",
                     column(width = 2, box(
                       checkboxInput("GCcluster", "Cluster", value = TRUE),
                       checkboxInput("theoreticalGC", "Normalize To Theoretical GC", value = FALSE),
                       htmlOutput("theoreticalGC"),
                       htmlOutput("GCspecies"), collapsible = TRUE, width = NULL, title = "Options"
                     ),
                     valueBoxOutput("GCboxP", width = NULL),
                     valueBoxOutput("GCboxW", width = NULL),
                     valueBoxOutput("GCboxF", width = NULL)),
                     box(h1("Per Sequence GC Content"),
                         h5("GC content (%) in sample, can either view total count or frequency"),
                         h5("Click sidebar on heatmap to change line plots"),
                         plotlyOutput("GCheatmap"),
                         br(),
                         plotlyOutput("GCSingle"),
                         width = 10)
             ),
             tabItem(tabName = "NC",
                     column(width = 2, box(
                       checkboxInput("Ncluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ), 
                     valueBoxOutput("NCboxP", width = NULL),
                     valueBoxOutput("NCboxW", width = NULL),
                     valueBoxOutput("NCboxF", width = NULL)),
                    box(
                      h1("Per base N content"),
                      h5("N content (%) in sample"),
                      h5("If dendrogram is truncated double click on dendrogram to resize"),
                      plotlyOutput("NCheatmap"),
                      br(),
                      plotlyOutput("NCsingle"),
                      width = 10)
             ),
             tabItem(tabName = "SLD", 
                     column(width = 2, box(
                       radioButtons(inputId="SLType", label="Value to plot",
                                    choices=c("Frequency","Counts"), selected = "Frequency"),
                       checkboxInput("SLcluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ),
                     valueBoxOutput("SLDboxP", width = NULL),
                     valueBoxOutput("SLDboxW", width = NULL),
                     valueBoxOutput("SLDboxF", width = NULL)),
                     box(
                       h1("Sequence Length Distribution"),
                       h5("Sequence length distribution in each sample, can either view total count or frequency"),
                       h5("Click sidebar on heatmap to change line plots"),
                       plotlyOutput("SLHeatmap"),
                       br(),
                       plotlyOutput("SLSingle"),
                       width = 10)
             ),
             tabItem(tabName = "SDL",
                     column(width = 2, box(
                       checkboxInput("Dupcluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ),
                     valueBoxOutput("SDLboxP", width = NULL),
                     valueBoxOutput("SDLboxW", width = NULL),
                     valueBoxOutput("SDLboxF", width = NULL)),
                     box(
                       h1("Sequence Duplication Levels"),
                       h5("Sequence duplication in each sample"),
                       h5("Click sidebar on heatmap to change line plots"),
                       plotlyOutput("DupHeatmap"),
                       br(),
                       plotlyOutput("DupSingle"), 
                       width = 10)
             ),
             tabItem(tabName = "OS",
                     column(width = 2, box(
                       checkboxInput("OScluster", "Cluster", value = TRUE), 
                       h5("Export Overrepresented Sequences"),
                       shinyDirButton(id = "dirOS", label = "Choose directory", title = ""),
                       collapsible = TRUE, width = NULL, title = "Options" 
                     ),
                     valueBoxOutput("OSboxP", width = NULL),
                     valueBoxOutput("OSboxW", width = NULL),
                     valueBoxOutput("OSboxF", width = NULL)),
                     box(
                       h1("Overrepresented Sequences"),
                       h5("Origin of Overrepresented sequences within each sample"),
                       plotlyOutput("OSummary"),
                       br(),
                       plotlyOutput("OSsingle"),
                       width = 10)
             ),
             tabItem(tabName = "AC",
                     column(width = 2, box(
                       selectInput("ACtype", "Choose Adapter Type",
                                   choices = c("Total",
                                               "Illumina Universal",
                                               "Illumina Small RNA",
                                               "Nextera Transposase")),
                       checkboxInput("ACcluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ),
                     valueBoxOutput("ACboxP", width = NULL),
                     valueBoxOutput("ACboxW", width = NULL),
                     valueBoxOutput("ACboxF", width = NULL)),
                     box(
                       h1("Adapter content"),
                       h5("Adapter content (%) across all reads"),
                       plotlyOutput("ACheatmap"),
                       br(),
                       plotlyOutput("ACsingle"),
                       width = 10)
             ),
             tabItem(tabName = "KC",
                     column(width = 2, box(
                       checkboxInput("KMcluster", "Cluster", value = TRUE), collapsible = TRUE, width = NULL, title = "Options"
                     ),
                     valueBoxOutput("KCboxP", width = NULL),
                     valueBoxOutput("KCboxW", width = NULL),
                     valueBoxOutput("KCboxF", width = NULL)),
                     box(
                       h1("Kmer Content"),
                       h5("Total Identified Kmer Count by Position.\nPlease select a file to see the top 6 Kmers."),
                       plotlyOutput("Kheatmap"),
                       br(),
                       plotlyOutput("Ksingle"),
                       width = 10)
             ),
             tabItem(tabName = "HTML",
                     column(width = 2, box(
                       radioButtons(inputId="omicsType", label="What type of -omic?",
                                    choices=c("Genome","Transcriptome"), selected = "Genome"),
                       htmlOutput("sequencedSpecies"),
                       h5("Output report for files"),
                       shinyDirButton(id = "dirs", label = "Choose directory", title = ""),
                       collapsible = TRUE, width = NULL, title = "Options"
                     )),
                     box(
                       h1("Output HTML Report Using the Default Template "),
                       h5("Select the type of data used in your study (-omic) and a closely related organism"),
                       h5("from the dropdown list. Upon selecting the applicable omic and species, select the directory"),
                       h5("containing the FASTQC files you wish to make the log for. Currently even if files have been loaded"),
                       h5("into the Shiny app the folder containing the data must still be selected."),
                       width = 10
                     )
             )
    )
    
  )
  
  ui <- dashboardPage(
    dashboardHeader(title = "ngsR::FASTQC"),
    dashboardSidebar(
      sidebarMenu(
        menuItem(text = "Summary", tabName = "BS"), 
        menuItem(text = "Total Sequences", tabName = "TS"),
        menuItemOutput("BQflag"),
        menuItemOutput("SQflag"),
        menuItemOutput("SCflag"),
        menuItemOutput("GCflag"),
        menuItemOutput("NCflag"),
        menuItemOutput("SLDflag"),
        menuItemOutput("SDLflag"),
        menuItemOutput("OSflag"),
        menuItemOutput("ACflag"),
        menuItemOutput("KCflag"),
        menuItem(text = "Output HTML Report", tabName = "HTML")
        
        
      ), width = 250
    ), body
  )
    
    
    server <- function(input, output, session){
      
      
    # set up reactives   
      values <- reactiveValues()
      
      #rective function repsonsible for loading in the selected files of just using the fdl supplied
      data <- reactive({
        volumes <- shinyFiles::getVolumes()
        shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session,
                                    filetypes = "zip")
        fileSelected <- shinyFiles::parseFilePaths(volumes, input$files)
        fileSelected <- as.character(fileSelected$datapath)
        selectedData <- getFastqcData(fileSelected)
        if(!length(input$files)){
          if(class(fastqcInput) != "FastqcDataList") selectedData <- getFastqcData(fastqcInput)
          else selectedData <- fastqcInput
        }
        selectedData
      })
      
      output$sequencedSpecies <- renderUI({
        
        if(input$omicsType == "Genome"){
          selectInput("omicSpecies", "Select species",
                      choices = genomes(gcTheoretical)$Name,
                      selected = "Hsapiens")
        }else{
          selectInput("omicSpecies", "Select species",
                      choices = transcriptomes(gcTheoretical)$Name,
                      selected = "Hsapiens")
        }
      })
      
      
      species <- observeEvent(input$omicSpecies, {
        values$omicSpecies <- input$omicSpecies
      })
      
      #export overrepresented
      
      expOS <- reactive({
        volumes <- shinyFiles::getVolumes()
        shinyFiles::shinyDirChoose(input, "dirOS", roots = volumes, session = session)
        dirSelected <- shinyFiles::parseDirPath(volumes, input$dirOS)
        as.character(dirSelected)
      })
      
       observe({
        if(length(expOS())){
        exportOverrepresented(data(), path = paste0(expOS(), "/OverrepSequences", "-", Sys.Date()), n = 10, noAdapters = TRUE)
        }
      })
      
   # export HTML    
      dir <- reactive({
        input$dirs
        volumes <- shinyFiles::getVolumes()
        shinyFiles::shinyDirChoose(input, "dirs", roots = volumes, session = session)
        dirSelected <- shinyFiles::parseDirPath(volumes, input$dirs)
        as.character(dirSelected)
      })
      
      observe({
        dir()
        if(length(dir())){
          withProgress(min = 0, max = 1, value = 0.8, message = "Writing report", {
            writeHtmlReport(dir(), species = values$omicSpecies, dataType = input$omicsType)
          })
          output$report2 <- renderText("Done!")
        }
      })
      
      
 # dynamic tabs to show pass warn Fail     
      
      output$BQflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          Category <- c()
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Per base sequence quality")
          
          items <- menuItemLogic(flags = flags)
        
          values$BQflag <- items[[1]]
          values$BQcolour <- items[[2]]
          values$BQcountF <- items[[3]]
          values$BQcountW <- items[[4]]
          values$BQcountP <- items[[5]]
          
          menuItem(text = "Base Quality", tabName = "BQ", badgeLabel = values$BQflag, badgeColor = values$BQcolour)
          
        }
        else{
          menuItem(text = "Base Quality", tabName = "BQ")
        }
      })
      
      
      output$SQflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Per sequence quality scores")
          
          items <- menuItemLogic(flags = flags)
          
          values$SQflag <- items[[1]]
          values$SQcolour <- items[[2]]
          values$SQcountF <- items[[3]]
          values$SQcountW <- items[[4]]
          values$SQcountP <- items[[5]]

          menuItem(text = "Sequence Quality", tabName = "SQ", badgeLabel = values$SQflag, badgeColor = values$SQcolour)
          
        }
        else{
          menuItem(text = "Sequence Quality", tabName = "SQ")
        }
      })
      
      output$SCflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Per base sequence content")
          
          items <- menuItemLogic(flags = flags)
          
          values$SCflag <- items[[1]]
          values$SCcolour <- items[[2]]
          values$SCcountF <- items[[3]]
          values$SCcountW <- items[[4]]
          values$SCcountP <- items[[5]]
          
          menuItem(text = "Sequence Content", tabName = "SC", badgeLabel = values$SCflag, badgeColor = values$SCcolour)
          
        }
        else{
          menuItem(text = "Sequence Content", tabName = "SC")
        }
      })
      
      output$GCflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Per sequence GC content")
          
          items <- menuItemLogic(flags = flags)
          
          values$GCflag <- items[[1]]
          values$GCcolour <- items[[2]]
          values$GCcountF <- items[[3]]
          values$GCcountW <- items[[4]]
          values$GCcountP <- items[[5]]
          
          menuItem(text = "GC Content", tabName = "GC", badgeLabel = values$GCflag, badgeColor = values$GCcolour)
          
        }
        else{
          menuItem(text = "GC Content", tabName = "GC")
        }
      })
      
      output$NCflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Per base N content")
          
          items <- menuItemLogic(flags = flags)
          
          values$NCflag <- items[[1]]
          values$NCcolour <- items[[2]]
          values$NCcountF <- items[[3]]
          values$NCcountW <- items[[4]]
          values$NCcountP <- items[[5]]
          
          menuItem(text = "N Content", tabName = "NC", badgeLabel = values$NCflag, badgeColor = values$NCcolour)
          
        }
        else{
          menuItem(text = "N Content", tabName = "NC")
        }
      })
      
      output$SLDflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Sequence Length Distribution")
          
          items <- menuItemLogic(flags = flags)
          
          values$SLDflag <- items[[1]]
          values$SLDcolour <- items[[2]]
          values$SLDcountF <- items[[3]]
          values$SLDcountW <- items[[4]]
          values$SLDcountP <- items[[5]]
          
          menuItem(text = "Sequence Length Distribution", tabName = "SLD", badgeLabel = values$SLDflag, badgeColor = values$SLDcolour)
          
        }
        else{
          menuItem(text = "Sequence Length Distribution", tabName = "SLD")
        }
      })
      
      output$SDLflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Sequence Duplication Levels")
          
          items <- menuItemLogic(flags = flags)
          
          values$SDLflag <- items[[1]]
          values$SDLcolour <- items[[2]]
          values$SDLcountF <- items[[3]]
          values$SDLcountW <- items[[4]]
          values$SDLcountP <- items[[5]]
          
          menuItem(text = "Sequence Duplicaiton Levels", tabName = "SDL", badgeLabel = values$SDLflag, badgeColor = values$SDLcolour)
          
        }
        else{
          menuItem(text = "Sequence Duplicaiton Levels", tabName = "SDL")
        }
      })
      
      output$OSflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Overrepresented sequences")
          
          items <- menuItemLogic(flags = flags)
          
          values$OSflag <- items[[1]]
          values$OScolour <- items[[2]]
          values$OScountF <- items[[3]]
          values$OScountW <- items[[4]]
          values$OScountP <- items[[5]]
          
          menuItem(text = "Overrepresented Sequences", tabName = "OS", badgeLabel = values$OSflag, badgeColor = values$OScolour)
          
        }
        else{
          menuItem(text = "Overrepresented Sequences", tabName = "OS")
        }
      })
      
      
      output$ACflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Adapter Content")
          
          items <- menuItemLogic(flags = flags)
          
          values$ACflag <- items[[1]]
          values$ACcolour <- items[[2]]
          values$ACcountF <- items[[3]]
          values$ACcountW <- items[[4]]
          values$ACcountP <- items[[5]]
          
          menuItem(text = "Adapter Content", tabName = "AC", badgeLabel = values$ACflag, badgeColor = values$ACcolour)
          
        }
        else{
          menuItem(text = "Adapter Content", tabName = "AC")
        }
      })
      
      output$KCflag <- renderMenu({
        if(!is.null(fastqcInput) | !is.null(input$files)){
          flags <- getSummary(data())
          flags <- subset(flags, flags$Category == "Kmer Content")
          
          items <- menuItemLogic(flags = flags)
          
          values$KCflag <- items[[1]]
          values$KCcolour <- items[[2]]
          values$KCcountF <- items[[3]]
          values$KCcountW <- items[[4]]
          values$KCcountP <- items[[5]]
          
          menuItem(text = "K-mer Content", tabName = "KC", badgeLabel = values$KCflag, badgeColor = values$KCcolour)
          
        }
        else{
          menuItem(text = "K-mer Content", tabName = "KC")
        }
      })
      
      #render fail value boxes
      output$BQboxF <- renderValBox(count = values$BQcountF, status = "FAIL", ic = "times", c = "red")
      
      output$SQboxF <- renderValBox(count = values$SQcountF, status = "FAIL", ic = "times", c = "red")
      
      output$SCboxF <- renderValBox(count = values$SCcountF, status = "FAIL", ic = "times", c = "red")
      
      output$GCboxF <- renderValBox(count = values$GCcountF, status = "FAIL", ic = "times", c = "red")
      
      output$NCboxF <- renderValBox(count = values$NCcountF, status = "FAIL", ic = "times", c = "red")
      
      output$SLDboxF <- renderValBox(count = values$SLDcountF, status = "FAIL", ic = "times", c = "red")
      
      output$SDLboxF <- renderValBox(count = values$SDLcountF, status = "FAIL", ic = "times", c = "red")
      
      output$OSboxF <- renderValBox(count = values$OScountF, status = "FAIL", ic = "times", c = "red")
      
      output$ACboxF <- renderValBox(count = values$ACcountF, status = "FAIL", ic = "times", c = "red")
      
      output$KCboxF <- renderValBox(count = values$KCcountF, status = "FAIL", ic = "times", c = "red")
      
      #render warn value boxes
      
      output$BQboxW <- renderValBox(count = values$BQcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$SQboxW <- renderValBox(count = values$SQcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$SCboxW <- renderValBox(count = values$SCcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$GCboxW <- renderValBox(count = values$GCcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$NCboxW <- renderValBox(count = values$NCcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$SLDboxW <- renderValBox(count = values$SLDcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$SDLboxW <- renderValBox(count = values$SDLcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$OSboxW <- renderValBox(count = values$OScountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$ACboxW <- renderValBox(count = values$ACcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      output$KCboxW <- renderValBox(count = values$KCcountW, status = "WARN", ic = "exclamation", c = "yellow")
      
      #render Pass value boxes
      
      output$BQboxP <- renderValBox(count = values$BQcountP, status = "PASS", ic = "check", c = "green")
      
      output$SQboxP <- renderValBox(count = values$SQcountP, status = "PASS", ic = "check", c = "green")
      
      output$SCboxP <- renderValBox(count = values$SCcountP, status = "PASS", ic = "check", c = "green")
      
      output$GCboxP <- renderValBox(count = values$GCcountP, status = "PASS", ic = "check", c = "green")
      
      output$NCboxP <- renderValBox(count = values$NCcountP, status = "PASS", ic = "check", c = "green")
      
      output$SLDboxP <- renderValBox(count = values$SLDcountP, status = "PASS", ic = "check", c = "green")
      
      output$SDLboxP <- renderValBox(count = values$SDLcountP, status = "PASS", ic = "check", c = "green")
      
      output$OSboxP <- renderValBox(count = values$OScountP, status = "PASS", ic = "check", c = "green")
      
      output$ACboxP <- renderValBox(count = values$ACcountP, status = "PASS", ic = "check", c = "green")
      
      output$KCboxP <- renderValBox(count = values$KCcountP, status = "PASS", ic = "check", c = "green")
      
      
      

      
      
# render plots       
      
      ####################
      # Summary
      ####################
      
      output$SummaryFlags <- renderPlotly({
        plotSummary(data(), usePlotly = TRUE,
                    cluster = input$Sumcluster, dendrogram = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      ####################
      # Read Totals
      ####################
      
      output$ReadTotals <- renderPlotly({
        plotReadTotals(data(), usePlotly = TRUE, duplicated = input$showDup) %>%
          layout(margin = list(l = 100, r = 200))
      })

      ####################
      # Base Quality
      ####################
            
      output$baseQualHeatmap <- renderPlotly({
        plotBaseQualities(data(),
                          usePlotly = TRUE,
                          plotType = "heatmap",
                          plotValue = input$BQplotValue,
                          cluster = input$BQcluster,
                          dendrogram = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      output$BaseQualitiesSingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        plotBaseQualities(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200, b = 50))
      })
      
      ####################
      # Sequence Quality
      ####################
            
      output$seqQualHeatmap <- renderPlotly({
        plotSequenceQualities(data(),
                              cluster = input$SQcluster,
                              counts = FALSE,
                              dendrogram = TRUE,
                              usePlotly = TRUE) %>% layout(margin = list(r = 200))
      })
      
      output$SeqQualitiesSingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        qualPlot <- plotSequenceQualities(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200),
                 legend = list(orientation = 'h', title = ""))
      })
      
      ####################
      # Sequence Content
      ####################   
      
      output$SCHeatmap <- renderPlotly({
        
        plotSequenceContent(data(),
                            cluster = input$SCcluster,
                            dendrogram = TRUE,
                            usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      
      output$SCsingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        plotSequenceContent(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      
      ####################
      # GC Content
      ####################
      
      output$theoreticalGC <- renderUI({
        if(input$theoreticalGC){
          radioButtons(inputId="theoreticalType", label="What type of data?",
                       choices=c("Genome","Transcriptome"), selected = "Genome")
        }
      })
      
      output$GCspecies <- renderUI({
        if(!is.null(input$theoreticalGC)){
          if(input$theoreticalGC){
            if(input$theoreticalType == "Genome"){
              selectInput("GCspecies", "Select species",
                          choices = genomes(gcTheoretical)$Name,
                          selected = "Hsapiens")
            }else{
              selectInput("GCspecies", "Select species",
                          choices = transcriptomes(gcTheoretical)$Name,
                          selected = "Hsapiens")
            }
          }
        }})
      
      
      
      output$GCheatmap <- renderPlotly({
        if(is.null(input$GCspecies)){
          GCspecies <- FALSE
        }else GCspecies <- input$GCspecies
        
        GCtype <- input$GCheatType == "Count"
        
        if(is.null(input$theoreticalGC)){
          plotGcContent(data(),
                        cluster = input$GCcluster,
                        plotType = "heatmap",
                        theoreticalType = input$theoreticalType,
                        dendrogram = TRUE,
                        usePlotly = TRUE) %>%
            layout(margin = list(r = 200))
        }else{
          plotGcContent(data(),
                        cluster = input$GCcluster,
                        plotType = "heatmap",
                        theoreticalType = input$theoreticalType,
                        theoreticalGC = input$theoreticalGC,
                        dendrogram = TRUE,
                        species = GCspecies,
                        usePlotly = TRUE) %>%
            layout(margin = list(r = 200))
        }
        
      })
      
      output$GCSingle <- renderPlotly({
        if(is.null(input$GCspecies)){
          GCspecies <- FALSE
        }else GCspecies <- input$GCspecies
        
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        if(is.null(input$theoreticalGC)){
          GCSingle <- plotGcContent(sub_fdl, usePlotly = TRUE,
                                    counts = FALSE)
        }else{
          GCSingle <- plotGcContent(sub_fdl, usePlotly = TRUE,
                                    counts = FALSE, theoreticalGC = input$theoreticalGC,
                                    theoreticalType = input$theoreticalType,
                                    species = GCspecies)
        }
        GCSingle %>%
          layout(margin = list(r = 200),
                 legend = list(orientation = 'h', title = ""))
      })
      
      
      ####################
      # N Content
      ####################
      
      output$NCheatmap <- renderPlotly({
        
        plotNContent(data(),
                     cluster = input$Ncluster,
                     dendrogram = TRUE,
                     usePlotly = TRUE) %>% layout(margin = list(r = 200))
      })
      
      
      # N Content single plot
      
      output$NCsingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        plotNContent(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200, b = 50))
      })
      
      
      ####################
      # Sequence Length Distribution
      ####################
      
      output$SLHeatmap <- renderPlotly({
        plotSequenceLengthDistribution(data(),
                                       cluster = input$SLcluster,
                                       dendrogram = TRUE,
                                       counts = input$SLType == "Counts",
                                       usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      output$SLSingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        plotSequenceLengthDistribution(sub_fdl, usePlotly = TRUE, plotType = "line") %>%
          layout(margin = list(r = 200))
      })
      
      ####################
      # Sequence Duplicaiton Levels
      ####################

      output$DupHeatmap <- renderPlotly({
        
        plotDuplicationLevels(data(),
                              cluster = input$Dupcluster,
                              dendrogram = TRUE,
                              usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      output$DupSingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        plotDuplicationLevels(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      ####################
      # Overrepresented sequences
      ####################
      
      output$OSummary <- renderPlotly({
        plotOverrepresentedSummary(data(),
                                   usePlotly = TRUE,
                                   cluster = input$OScluster,
                                   dendrogram = TRUE) %>%
          layout(margin = list(r = 200))
        
      })
      
      
      output$OSsingle <- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        plotOverrepresentedSummary(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      })
      
      ####################
      # Adapter Content
      ####################     
     
      output$ACheatmap <- renderPlotly({
        ACplot <- plotAdapterContent(data(),
                                     adapterType = input$ACtype,
                                     usePlotly = TRUE,
                                     dendrogram = TRUE,
                                     cluster = input$ACcluster)
        if(!is.null(ACplot)) ACplot %>% layout(margin = list(r = 200))
        else stop(paste("Sequences did not contain any", input$ACtype, "content, please select another."))
      })
      
      output$ACsingle<- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(fileName(data()) == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        ACsing <- plotAdapterContent(sub_fdl, usePlotly = TRUE)
        
        if(!is.null(ACsing)) ACsing %>%
          layout(margin = list(r = 200),
                 legend = list(orientation = 'h', title = ""))
        else stop(paste("Sequences did not contain any",
                        input$ACtype, "content, please select another."))
      })
      
      
      ####################
      # k-mer Content
      ####################
       
      output$Kheatmap <- renderPlotly({
        Kplot <- plotKmers(data(),
                           usePlotly = TRUE,
                           cluster = input$KMcluster,
                           dendrogram = TRUE)
        if(!is.null(Kplot)) Kplot %>% layout(margin = list(r = 200))
        else stop(paste("Samples have no Kmer content"))
      })
      
      output$Ksingle<- renderPlotly({
        if(is.null(event_data("plotly_click")$key[[1]])){
          num <- 1
        }else {
          click <- event_data("plotly_click")
          num <- which(grepl(click$key[[1]], fileName(data())))
        }
        sub_fdl <- data()[[num]]
        Ksing <- plotKmers(sub_fdl, usePlotly = TRUE)
        
        if(!is.null(Ksing)) Ksing %>%
          layout(margin = list(r = 200, b = 50))
        else stop(paste("Library did not contain any identified Kmers please select another."))
      })
    }

    shinyApp(ui = ui, server = server)

}
