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
#' @importFrom reshape2 dcast
#' @importFrom magrittr %>%
#' @importFrom plotly layout
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom plotly event_data
#' @importFrom shiny shinyUI
#' @importFrom shiny fluidPage
#' @importFrom shiny navbarPage
#' @importFrom shiny tabPanel
#' @importFrom shiny splitLayout
#' @importFrom shiny fixedPanel
#' @importFrom shiny absolutePanel
#' @importFrom shiny h1
#' @importFrom shiny h5
#' @importFrom shiny observe
#' @importFrom shiny sidebarPanel
#' @importFrom shiny radioButtons
#' @importFrom shiny checkboxInput
#' @importFrom shiny htmlOutput
#' @importFrom shiny selectInput
#' @importFrom shiny plotOutput
#' @importFrom shiny sliderInput
#' @importFrom shiny renderUI
#' @importFrom shiny renderPrint
#' @importFrom shiny shinyApp
#' @importFrom shiny textOutput
#' @importFrom shiny renderText
#' @importFrom shiny reactive
#' @importFrom shiny withProgress
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles shinyDirChoose
#' @importFrom shinyFiles shinyDirButton
#' @importFrom shinyFiles parseFilePaths
#' @importFrom shinyFiles parseDirPath
#' @importFrom shinyFiles getVolumes
#'
#' @export
#' @rdname fastqcShiny
fastqcShiny <- function(fastqcInput = NULL){

  if(!is.null(fastqcInput)){
    stopifnot(length(fastqcInput) > 1)
    stopifnot(class(fastqcInput) == "character" | class(fastqcInput)[1] == "FastqcDataList")
  }


  ui <- shinyUI(
    fluidPage(
      navbarPage(
        "ngsReports::FASTQC",
        #first panel is summary
        tabPanel(
          "fastQC Flags Summary",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                h5("Choose Fastqc Report:"),
                shinyFiles::shinyFilesButton(id = "files", label = "Choose files", multiple = TRUE, title = ""),
                h5(""),
                textOutput("report"),
                h5(""),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Summary of fastQC Flags"),
              h5("Heatmap of fastQC flags (pass, warning or fail) for each fastQC report"),
              plotlyOutput("SummaryFlags"),
              width = "70%", left = "30%", right = "0%", height = "1100px")
          )
        ),
        tabPanel(
          "Total Sequences",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Read Totals"),
              h5("Total number of unique and duplicated reads in each sample"),
              plotlyOutput("ReadTotals"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Duplicated sequences",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                checkboxInput("Dupcluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Sequence Duplication levels"),
              h5("Sequence duplication in each sample"),
              h5("Click sidebar on heatmap to change line plots"),
              plotlyOutput("DupHeatmap"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Per Base Sequence Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Per Base Sequence Content"),
              h5("Per base sequence content in each sample, colours at each base indicate sequence bias"),
              h5("1 - G = opacity, A = Green, T = Red, C = Blue"),
              h5("if each base is equally represented then should be light grey"),
              plotlyOutput("SCHeatmap"),
              plotlyOutput("SCsingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px"))),
        tabPanel(
          "Per Base Sequence Quality",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="BQplotValue", label="Base Quality",
                             choices=c("Mean","Median"), selected = "Mean"),
                checkboxInput("BQcluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Base Quality"),
              h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
              h5("Click sidebar on heatmap to change line plots"),
              plotlyOutput("baseQualHeatmap"),
              plotlyOutput("BaseQualitiesSingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px"))),
        tabPanel(
          "Per Sequence Quality Scores",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="SQType", label="Sequence Quality",
                             choices=c("Frequency","Counts"), selected = "Frequency"),
                checkboxInput("SQcluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Sequence Quality"),
              h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
              h5("Click sidebar on heatmap to change line plots"),
              plotlyOutput("seqQualHeatmap"),
              plotlyOutput("SeqQualitiesSingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px"))),
        tabPanel(
          "% GC Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="theoreticalType", label="What type of data?",
                             choices=c("Genome","Transcriptome"), selected = "Genome"),
                radioButtons(inputId="GCheatType", label="Value to plot",
                             choices=c("Frequency","Count"), selected = "Frequency"),
                checkboxInput("GCcluster", "Cluster Filenames", value = TRUE),
                htmlOutput("theoreticalGC"),
                htmlOutput("GCspecies"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("GC content in reads"),
              h5("GC content (%) in sample, can either view total count or frequency"),
              h5("Click sidebar on heatmap to change line plots"),
              plotlyOutput("GCheatmap"),
              plotlyOutput("GCSingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px"))),
        tabPanel(
          "Sequence Length Distribution",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="SLType", label="Value to plot",
                             choices=c("Frequency","Counts"), selected = "Frequency"),
                checkboxInput("SLcluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Sequence Length Distribution for all reads"),
              h5("Sequence length distribution in each sample, can either view total count or frequency"),
              h5("Click sidebar on heatmap to change line plots"),
              plotlyOutput("SLHeatmap"),
              plotlyOutput("SLSingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px"))),
        tabPanel(
          "Overrepresented Sequences",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                checkboxInput("ORcluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Overrepresented Sequences"),
              h5("Origin of Overrepresented sequences within each sample"),
              plotlyOutput("OSummary"),
              width = "70%", left = "30%", right = "0%", height = "1100px"))),
        tabPanel(
          "% N-Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                checkboxInput("Ncluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("N content in reads"),
              h5("N content (%) in sample"),
              h5("If dendrogram is truncated double click on dendrogram to resize"),
              plotlyOutput("NCheatmap"),
              width = "70%", left = "30%", right = "0%", height = "1100px")
          )
        ),
        tabPanel(
          "Adapter Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                selectInput("ACtype", "Choose Adapter Type",
                            choices = c("Total",
                                        "Illumina Universal",
                                        "Illumina Small RNA",
                                        "Nextera Transposase")),
              checkboxInput("ACcluster", "Cluster Filenames", value = TRUE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Adapter content"),
              h5("Adapter content (%) across all reads"),
              plotlyOutput("ACheatmap"),
              plotlyOutput("ACsingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px")
          )
        ),
        tabPanel(
          "Kmer Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Kmer Content"),
              h5("Total Identified Kmer Count by Position.\nPlease select a file to see the top 6 Kmers."),
              plotlyOutput("Kheatmap"),
              plotlyOutput("Ksingle"),
              width = "70%", left = "30%", right = "0%", height = "1100px")
          )
        ),
        tabPanel(
          "Output HTML Report",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="omicsType", label="What type of -omic?",
                             choices=c("Genome","Transcriptome"), selected = "Genome"),
                htmlOutput("sequencedSpecies"),
                h5("Output report for files"),
                shinyDirButton(id = "dirs", label = "Choose directory", title = ""),
                textOutput("report2"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Output HTML Report Using the Default Template "),
              h5("Select the type of data used in your study (-omic) and a closely related organism"),
              h5("from the dropdown list. Upon selecting the applicable omic and species, select the directory"),
              h5("containing the FASTQC files you wish to make the log for. Currently even if files have been loaded"),
              h5("into the Shiny app the folder containing the data must still be selected."),
              width = "70%", left = "30%", right = "0%")
      )
    )
      )))

  server <- function(input, output, session){


    #rective function repsonsible for loading in the selected files of just using the fdl supplied
    data <- reactive({
      volumes <- shinyFiles::getVolumes()
      shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session,
                                  filetypes = "zip")
      fileSelected <- shinyFiles::parseFilePaths(volumes, input$files)
      fileSelected <- as.character(fileSelected$datapath)
      selectedData <- getFastqcData(fileSelected)
      if(is.null(input$files)){
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


    species <- reactive({
      input$omicSpecies
    })

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
          writeHtmlReport(dir(), species = species(), dataType = input$omicsType)
        })
        output$report2 <- renderText("Done!")
      }
    })





# Summary heatmap in first tab
    output$SummaryFlags <- renderPlotly({
      plotSummary(data(), usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    })


# report the number of files you ve read in
    output$report <- renderText({
      if(!length(data())){
        print("No Files Loaded")
      }else{
        paste(length(data()), "Files Loaded")
      }
    })


    # plot read totals

    output$ReadTotals <- renderPlotly({
      # plotDeduplicatedTotalsPlotly(data(), bars = input$RDLbar) %>%
      #   layout(margin = list(r = 200))
      plotReadTotals(data(), usePlotly = TRUE, duplicated = TRUE) %>%
        layout(margin = list(r = 200))
    })

    output$OSummary <- renderPlotly({
      plotOverrepresentedSummary(data(),
                                 usePlotly = TRUE) %>%
        layout(margin = list(r = 200))

    })

    output$DupHeatmap <- renderPlotly({

      plotDuplicationLevels(data(),
                            clusterNames = input$Dupcluster,
                            dendrogram = TRUE,
                            usePlotly = TRUE) %>%
        layout(margin = list(r = 200, l = 0))
    })

    output$SCHeatmap <- renderPlotly({

      plotSequenceContent(data(),
                            usePlotly = TRUE) %>%
        layout(margin = list(r = 200, l = 0))
    })


    output$SCsingle <- renderPlotly({
      if(is.null(event_data("plotly_click")$key[[1]])){
        num <- 1
      }else {
        click <- event_data("plotly_click")
        num <- which(fileName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotSequenceContent(sub_fdl, usePlotly = TRUE, plotType = "line") %>%
        layout(margin = list(r = 200, l = 0))
    })


    # start rendering the input buttons for
    output$theoreticalGC <- renderUI({
      if(input$GCheatType == "Frequency"){
        checkboxInput("theoreticalGC", "Normalize Using Theoretical GC", value = FALSE)
      }
    })


    output$GCspecies <- renderUI({
      if(!is.null(input$theoreticalGC)){
      if(input$theoreticalGC){
        if(input$theoreticalType == "Genome"){
        selectInput("GCspecies", "Select species for Theoretical GC",
                    choices = genomes(gcTheoretical)$Name,
                    selected = "Hsapiens")
      }else{
        selectInput("GCspecies", "Select species for Theoretical GC",
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
                      clusterNames = input$GCcluster,
                      plotType = "heatmap",
                      theoreticalType = input$theoreticalType,
                      dendrogram = TRUE,
                      usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }else{
        plotGcContent(data(),
                      clusterNames = input$GCcluster,
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
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
    })


    output$SLHeatmap <- renderPlotly({
      plotSequenceLengthDistribution(data(),
                                     clusterNames = input$SLcluster,
                                     dendrogram = TRUE,
                                     counts= FALSE,
                                     usePlotly = TRUE) %>%
        layout(margin = list(r = 200, l = 0))
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
        layout(margin = list(r = 200, l = 100))
    })

    output$baseQualHeatmap <- renderPlotly({
      plotBaseQualities(data(),
                        usePlotly = TRUE,
                        plotType = "heatmap",
                        plotValue = input$BQplotValue,
                        clusterNames = input$BQcluster,
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
        layout(margin = list(r = 200, l = 100))
    })


    output$seqQualHeatmap <- renderPlotly({
      plotSequenceQuality(data(),
                                   clusterNames = input$SQcluster,
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
      sub_fdl <- data()[num]
      qualPlot <- plotSequenceQuality(sub_fdl, usePlotly = TRUE, plotType = "line") %>%
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
    })

    output$NCheatmap <- renderPlotly({

        plotNContentPlotly(data(),
                           clusterNames = input$Ncluster,
                           dendrogram = TRUE,
                           usePlotly = TRUE) %>% layout(margin = list(r = 200))
      })

    # By default, this method will now plot the sum of all adapter types
    # Do we need to change the drop-down menu to accept the option "all"?
    output$ACheatmap <- renderPlotly({
      ACplot <- plotAdapterContent(data(),
                         adapterType = input$ACtype,
                         usePlotly = TRUE,
                         clusterNames = input$ACcluster)
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
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
      else stop(paste("Sequences did not contain any",
                      input$ACtype, "content, please select another."))
    })

    output$Kheatmap <- renderPlotly({
      Kplot <- plotKmers(data(),
                                   usePlotly = TRUE)
      if(!is.null(Kplot)) Kplot %>% layout(margin = list(r = 200))
      else stop(paste("Samples have no Kmer content"))
    })

    output$Ksingle<- renderPlotly({
      if(is.null(event_data("plotly_click")$key[[1]])){
        num <- 1
      }else {
        click <- event_data("plotly_click")
        num <- which(fileName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      Ksing <- plotKmers(sub_fdl, usePlotly = TRUE)

      if(!is.null(Ksing)) Ksing %>%
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
      else stop(paste("Library did not contain any identified Kmers please select another."))
    })


  }

  shinyApp(ui = ui, server = server)

}

