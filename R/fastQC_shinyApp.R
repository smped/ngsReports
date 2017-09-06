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
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles parseFilePaths
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
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Summary of fastQC Flags"),
              h5("Heatmap of fastQC flags (pass, warning or fail) for each fastQC report"),
              plotlyOutput("SummaryFlags"),
              width = "70%", left = "30%", right = "0%")
          )
        ),
        tabPanel(
          "Total Reads",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="RDLbar", label="Bar presentation",
                             choices=c("stacked","adjacent"), selected = "stacked"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Read Totals"),
              h5("Total number of unique and duplicated reads in each sample"),
              plotlyOutput("ReadTotals"),
              # plotlyOutput("RDSingle"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Per Base Sequence Quality",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="BQplotValue", label="Base Quality",
                             choices=c("Mean","Median"), selected = "Mean"),
                checkboxInput("BQcluster", "Cluster Filenames", value = FALSE),
                htmlOutput("BQdendro"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Base Quality"),
              h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
              h5("Click sidebar on heatmap to change line plots"),
              h5("If dendrogram is truncated double click on dendrogram to resize"),
              plotlyOutput("baseQualHeatmap"),
              plotlyOutput("BaseQualitiesSingle"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Per Sequence Quality Scores",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="SQType", label="Sequence Quality",
                             choices=c("Frequency","Counts"), selected = "Frequency"),
                checkboxInput("SQcluster", "Cluster Filenames", value = FALSE),
                htmlOutput("SQdendro"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Sequence Quality"),
              h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
              h5("Click sidebar on heatmap to change line plots"),
              h5("If dendrogram is truncated double click on dendrogram to resize"),
              plotlyOutput("seqQualHeatmap"),
              plotlyOutput("SeqQualitiesSingle"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "% GC Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="theoreticalType", label="What type of data?",
                             choices=c("Genome","Transcriptome"), selected = "Genome"),
                radioButtons(inputId="GCheatType", label="Value to plot",
                             choices=c("Frequency","Count"), selected = "Frequency"),
                checkboxInput("GCcluster", "Cluster Filenames", value = FALSE),
                htmlOutput("GCdendro"),
                htmlOutput("theoreticalGC"),
                htmlOutput("GCspecies"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("GC content in reads"),
              h5("GC content (%) in sample, can either view total count or frequency"),
              h5("Click sidebar on heatmap to change line plots"),
              h5("If dendrogram is truncated double click on dendrogram to resize"),
              plotlyOutput("GCheatmap"),
              plotlyOutput("GCSingle"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Sequence Length Distribution",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="SLType", label="Value to plot",
                             choices=c("Frequency","Count"), selected = "Frequency"),
                checkboxInput("SLcluster", "Cluster Filenames", value = FALSE),
                htmlOutput("SLdendro"),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Sequence Length Distribution for all reads"),
              h5("Sequence length distribution in each sample, can either view total count or frequency"),
              h5("Click sidebar on heatmap to change line plots"),
              h5("If dendrogram is truncated double click on dendrogram to resize"),
              plotlyOutput("SLHeatmap"),
              plotlyOutput("SLSingle"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Overrepresented Sequences",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="ORType", label="Individual or Overall",
                             choices=c("Individual","Overall"), selected = "Overall"),
                checkboxInput("ORcluster", "Cluster Filenames", value = FALSE),
                sliderInput("ORslide", "Number of seq", min = 1, max = 20, value = 10),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Overrepresented Sequences"),
              h5("Overrepresented sequences in each sample, can either view sequence on an individual or overall basis"),
              plotlyOutput("overRepHeatmap"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "% N-Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                checkboxInput("Ncluster", "Cluster Filenames", value = FALSE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("N content in reads"),
              h5("N content (%) in sample"),
              h5("If dendrogram is truncated double click on dendrogram to resize"),
              plotlyOutput("NCheatmap"),
              width = "70%", left = "30%", right = "0%")
          )
        ),
        tabPanel(
          "Adapter Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                selectInput("ACtype", "Choose Adapter Type",
                            choices = c("Illumina Universal",
                                        "Illumina Small RNA",
                                        "Nextera Transposase")),
              # checkboxInput("ACcluster", "Cluster Filenames", value = FALSE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Adapter content"),
              h5("Adapter content (%) across all reads"),
              plotlyOutput("ACheatmap"),
              plotlyOutput("ACsingle"),
              width = "70%", left = "30%", right = "0%")
          )
        ),
        tabPanel(
          "Kmer Content",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                # checkboxInput("ACcluster", "Cluster Filenames", value = FALSE),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Kmer Content"),
              h5("-log(10) P-value for Kmers"),
              plotlyOutput("Kheatmap"),
              plotlyOutput("Ksingle"),
              width = "70%", left = "30%", right = "0%")
          )
        )
      )
    )
  )

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

# Summary heatmap in first tab
    output$SummaryFlags <- renderPlotly({
      plotSummary(fdl, usePlotly = TRUE) %>%
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
      plotReadTotals(data(), usePlotly = TRUE, duplicated = TRUE, bars = input$RDLbar) %>%
        layout(margin = list(r = 200))
    })

    # output$RDSingle <- renderPlotly({
    #   if(is.null(event_data("plotly_click")$key[[1]])){
    #     num <- 1
    #   }else {
    #     click <- event_data("plotly_click")
    #     num <- which(fileName(data()) == click$key[[1]])
    #   }
    #   sub_fdl <- data()[num]
    #   plotDuplicationLevels(sub_fdl, type = "Total sequences", usePlotly = TRUE) %>%
    #     layout(margin = list(r = 200),
    #            legend = list(orientation = 'h', title = ""))
    # })


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



    output$GCdendro <- renderUI({
      if(input$GCcluster) {
        checkboxInput("GCdendro", "Plot Dendrogram?", value = FALSE)
      }
    })



    output$GCheatmap <- renderPlotly({
      if(is.null(input$GCspecies)){
        GCspecies <- FALSE
      }else GCspecies <- input$GCspecies

      GCtype <- input$GCheatType == "Count"

      if(is.null(input$GCdendro)){
        GCdendro <- FALSE
      }else{
        GCdendro <- input$GCdendro
      }

      if(is.null(input$theoreticalGC)){
        plotGcContent(data(),
                      clusterNames = input$GCcluster,
                      plotType = "heatmap",
                      theoreticalType = input$theoreticalType,
                      dendrogram = GCdendro,
                      usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }else{
        plotGcContent(data(),
                      clusterNames = input$GCcluster,
                      plotType = "heatmap",
                      theoreticalType = input$theoreticalType,
                      theoreticalGC = input$theoreticalGC,
                      dendrogram = GCdendro,
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
        #I've broken this sorry. Not sure how to fix this
        num <- which(fileName(data()) == click$key[[1]])
      }
      GCtype <- input$GCheatType == "Count"
      sub_fdl <- data()[[num]]
      if(is.null(input$theoreticalGC)){
        GCSingle <- plotGcContent(sub_fdl, usePlotly = TRUE,
                                  counts = GCtype)
      }else{
        GCSingle <- plotGcContent(sub_fdl, usePlotly = TRUE,
                                  counts = GCtype, theoreticalGC = input$theoreticalGC,
                                  theoreticalType = input$theoreticalType,
                                  species = GCspecies)
      }
      GCSingle %>%
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
    })


    output$SLdendro <- renderUI({
      if(input$SLcluster) {
        checkboxInput("SLdendro", "Plot Dendrogram?", value = FALSE)
      }
    })

    output$SLHeatmap <- renderPlotly({
      if(is.null(input$SLdendro)){
        SLdendro <- FALSE
      }else{
        SLdendro <- input$SLdendro
      }
      SLtype <- input$SLType == "Counts"
      plotSequenceLengthDistribution(data(),
                                     clusterNames = input$SLcluster,
                                     dendrogram = SLdendro,
                                     counts= SLtype,
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

    output$overRepHeatmap <- renderPlotly({
      plotOverrepresentedHeatmapPlotly(data(),
                                       clusterNames = input$ORcluster,
                                       method = input$ORType,
                                       nSeq = input$ORslide,
                                       usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    })

    output$BQdendro <- renderUI({
      if(input$BQcluster) {
        checkboxInput("BQdendro", "Plot Dendrogram?", value = FALSE)
      }
    })

    output$baseQualHeatmap <- renderPlotly({
      if(is.null(input$BQdendro)){
        BQdendro <- FALSE
      }else{
        BQdendro <- input$BQdendro
      }
      plotBaseQualities(data(),
                        usePlotly = TRUE,
                        plotType = "heatmap",
                        plotValue = input$BQplotValue,
                        clusterNames = input$BQcluster,
                        dendrogram = BQdendro) %>%
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

    output$SQdendro <- renderUI({
      if(input$SQcluster) {
        checkboxInput("SQdendro", "Plot Dendrogram?", value = FALSE)
      }
    })

    output$seqQualHeatmap <- renderPlotly({
      SQtype <- input$SQType == "Counts"
      if(is.null(input$SQdendro)){
        SQdendro <- FALSE
      }else{
        SQdendro <- input$SQdendro
      }
      plotSequenceQualitiesHeatmap(data(),
                                   clusterNames = input$SQcluster,
                                   type = SQtype,
                                   dendrogram = SQdendro,
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
      qualPlot <- plotSequenceQualities(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
    })

    output$NCdendro <- renderUI({
      if(input$Ncluster) {
        checkboxInput("NCdendro", "Plot Dendrogram?", value = FALSE)
      }
    })

    output$NCheatmap <- renderPlotly({
      if(is.null(input$NCdendro)){
        NCdendro <- FALSE
      }else{
        NCdendro <- input$NCdendro
      }
        plotNContentPlotly(data(),
                           clusterNames = input$Ncluster,
                           dendrogram = NCdendro,
                           usePlotly = TRUE) %>% layout(margin = list(r = 200))
      })

    output$ACheatmap <- renderPlotly({
      ACplot <- plotAdapterContent(data(),
                         adapterType = input$ACtype,
                         usePlotly = TRUE)
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
      sub_fdl <- data()[num]
      ACsing <- plotAdapterContent(sub_fdl,
                                   adapterType = input$ACtype,
                                   plotType = "line",
                                   usePlotly = TRUE)

      if(!is.null(ACsing)) ACsing %>%
        layout(margin = list(r = 200, l = 100),
               legend = list(orientation = 'h', title = ""))
      else stop(paste("Sequences did not contain any",
                      input$ACtype, "content, please select another."))
    })

    output$Kheatmap <- renderPlotly({
      Kplot <- plotKmerHeatmap(data(),
                                   usePlotly = TRUE)
      if(!is.null(Kplot)) Kplot %>% layout(margin = list(r = 200))
      else stop(paste("Samples have no Kmer content"))
    })


  }

  shinyApp(ui = ui, server = server)

}

