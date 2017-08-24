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
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles parseFilePaths
#' @importFrom shinyFiles getVolumes
#'
#' @export
#' @rdname fastqcShiny
fastqcShiny <- function(fastqcInput = NULL){

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
              shinyFiles::shinyFilesButton(id = "files", label = "Choose files", multiple = TRUE, title = "Please select a file:"),
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
              h1("Duplication levels in reads"),
              h5("Total number of unique and duplicated reads in each sample"),
              plotlyOutput("ReadDuplication"),
              width = "70%", left = "30%", right = "0%"))),
        tabPanel(
          "Per Base Sequence Quality",
          splitLayout(
            fixedPanel(
              sidebarPanel(
                radioButtons(inputId="BQType", label="Base Quality",
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
                checkboxInput("GCtheory", "Normalize using theoretical GC", value = FALSE),
                htmlOutput("GCspecies"),
                radioButtons(inputId="GCheatType", label="Value to plot",
                             choices=c("Count","Frequency"), selected = "Frequency"),
                checkboxInput("GCcluster", "Cluster Filenames", value = FALSE),
                htmlOutput("GCdendro"),
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
        )
      )
    )
  )

  server <- function(input, output, session){


    #This function is repsonsible for loading in the selected file
    data <- reactive({
      volumes <- shinyFiles::getVolumes()
      shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session)
        fileSelected <- shinyFiles::parseFilePaths(volumes, input$files)
        fileSelected <- as.character(fileSelected$datapath)
        selectedData <- getFastqcData(fileSelected)
        if(is.null(input$files)){
          if(class(fastqcInput) != "FastqcDataList") selectedData <- getFastqcData(fastqcInput)
          else selectedData <- fastqcInput
        }
        selectedData
    })

    output$SummaryFlags <- renderPlotly({
      plotSummary(data(), usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    })


    output$ReadDuplication <- renderPlotly({
      plotDeduplicatedTotalsPlotly(data(), bars = input$RDLbar) %>%
        layout(margin = list(r = 200))
    })

    output$GCspecies <- renderUI({
      if(input$GCtheory){
        selectInput("GCspecies", "Select species for Theoretical GC",
                    choices = genomes(gcTheoretical)$Name,
                    selected = "Hsapiens")
      }
    })

    output$GCdendro <- renderUI({
      if(input$GCcluster) {
        checkboxInput("GCdendro", "Plot Dendrogram?", value = FALSE)
      }
    })

    output$GCheatmap <- renderPlotly({
      GCtype <- input$GCheatType == "Count"
      if(is.null(input$GCdendro)){
        plotGcHeatmap(data(),
                      clusterNames = input$GCcluster,
                      counts = GCtype,
                      GCtheory = input$GCtheory,
                      species = input$GCspecies,
                      usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }else{
        plotGcHeatmap(data(),
                      clusterNames = input$GCcluster,
                      counts = GCtype,
                      GCtheory = input$GCtheory,
                      dendrogram = input$GCdendro,
                      species = input$GCspecies,
                      usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }

    })

    output$GCSingle <- renderPlotly({
      if(is.null(event_data("plotly_click")$key[[1]])){
        num <- 1
      }else {
        click <- event_data("plotly_click")
        num <- which(fileName(data()) == click$key[[1]])
      }
      GCtype <- input$GCheatType == "Count"
      sub_fdl <- data()[num]
      plotGcContent(sub_fdl, usePlotly = TRUE, counts = GCtype) %>%
        layout(margin = list(r = 200, l = 100),
      legend = list(orientation = 'h', title = ""))
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
        plotBaseQualitiesHeatmap(data(),
                                 clusterNames = input$BQcluster,
                                 type = input$BQType,
                                 usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }else{
        plotBaseQualitiesHeatmap(data(),
                                 clusterNames = input$BQcluster,
                                 type = input$BQType,
                                 dendrogram = input$BQdendro,
                                 usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }

    })

    output$BaseQualitiesSingle <- renderPlotly({
      if(is.null(event_data("plotly_click")$key[[1]])){
        num <- 1
      }else {
        click <- event_data("plotly_click")
        num <- which(fileName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[num]
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
        plotSequenceQualitiesHeatmap(data(),
                                     clusterNames = input$SQcluster,
                                     type = SQtype,
                                     usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }else{
        plotSequenceQualitiesHeatmap(data(),
                                     clusterNames = input$SQcluster,
                                     type = SQtype,
                                     dendrogram = input$SQdendro,
                                     usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }
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
        plotNContentPlotly(data(),
                           clusterNames = input$Ncluster,
                           usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }else{
        plotNContentPlotly(data(),
                           clusterNames = input$Ncluster,
                           dendrogram = input$NCdendro,
                           usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }
    })

  }

  shinyApp(ui = ui, server = server)

}

