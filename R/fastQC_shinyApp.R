#' Run fastQC shiny app
#'
#' @description Returns a shiny app interface to parse many fastQC objects
#'
#' @details Currently some plots can take a while to render if the \code{FastqcDataList} passed to
#' \code{fastqcInput} has many elements
#'
#' @param fastqcInput can be a \code{FastqcFileList}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files.
#' @param subsetAll a \code{character} vector of length 1 to subset all files input into shiny app by
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
#'
#' @export
#' @rdname fastqcShiny
fastqcShiny <- function(fastqcInput, subsetAll = ""){

  # Shouldn't this just be getFastqcData(fastqcInput)?
  if(class(fastqcInput) != "FastqcDataList"){
    fdl <- getFastqcData(fastqcInput)
  }
  if(class(fastqcInput) == "FastqcDataList"){
    fdl <- fastqcInput
  }

  fdl <- fdl[grepl(subsetAll, fileNames(fdl))]

  ui <- shinyUI(
    fluidPage(
      navbarPage(
        "fastqcR",
        #first panel is summary
        tabPanel(
          "fastQC Flags Summary",
          splitLayout(
            fixedPanel(),
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
                selectInput("BQheight", "Plot Height", choices = c("auto", 250, 500, 1000)),
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("Base Quality"),
              h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
              plotlyOutput("baseQualHeatmap"),
              plotOutput("baseQualIndv"),
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
              h1("Base Quality"),
              h5("Per base sequence quality in each sample, can either view mean or median for each cycle"),
              plotlyOutput("seqQualHeatmap"),
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
                width = "20%", left = "0%", right = "80%"
              ), width = "20%"),
            absolutePanel(
              h1("GC content in reads"),
              h5("GC content (%) in sample, can either view total count or frequency"),
              plotlyOutput("GCheatmap"),
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
              plotlyOutput("NCheatmap"),
              width = "70%", left = "30%", right = "0%")
          )
        )
      )
    )
  )

  server <- function(input, output){
    output$SummaryFlags <- renderPlotly({
      plotSummary(fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    })


    output$ReadDuplication <- renderPlotly({
      plotDeduplicatedTotalsPlotly(fdl, bars = input$RDLbar) %>%
        layout(margin = list(r = 200))
    })

    output$GCspecies <- renderUI({
      if(input$GCtheory){
        selectInput("GCspecies", "Select species for Theoretical GC",
                    choices = genomes(gcTheoretical),
                    selected = "Hsapiens")
      }
    })


    output$GCheatmap <- renderPlotly({
      GCtype <- input$GCheatType == "Count"
      plotGCHeatmapPlotly(fdl,
                          clusterNames = input$GCcluster,
                          counts = GCtype,
                          GCtheory = input$GCtheory,
                          species = input$GCspecies,
                          usePlotly = TRUE) %>%
        layout(margin = list(r = 200))

    })

    output$overRepHeatmap <- renderPlotly({
      plotOverrepresentedHeatmapPlotly(fdl,
                                       clusterNames = input$ORcluster,
                                       method = input$ORType,
                                       nSeq = input$ORslide) %>%
        layout(margin = list(r = 200))
    })

    output$BQdendro <- renderUI({
      if(input$BQcluster) {
        checkboxInput("BQdendro", "plot Dendrogram", value = FALSE)
      }
    })

    output$baseQualHeatmap <- renderPlotly({
      if(is.null(input$BQdendro)){
        plotBaseQualitiesPlotly(fdl,
                                clusterNames = input$BQcluster,
                                type = input$BQType,
                                usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }else{
        plotBaseQualitiesPlotly(fdl,
                                clusterNames = input$BQcluster,
                                type = input$BQType,
                                dendrogram = input$BQdendro,
                                usePlotly = TRUE) %>% layout(margin = list(r = 200))
      }

    })

    output$click <- renderPrint({d <- event_data("plotly_click")
    d$key[[1]]
    })



    output$baseQualIndv <- shiny::renderPlot({
      click <- event_data("plotly_click")
      key <- click$key[[1]]

      if(!is.null()){plotBaseQualities(FastqcFile())}
    })


    output$seqQualHeatmap <- renderPlotly({
      SQtype <- input$SQType == "Counts"
      plotSequenceQualitiesHeatmap(fdl,
                                   counts = SQtype,
                                   clusterNames = input$SQcluster,
                                   usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    })

    output$NCheatmap <- renderPlotly({
      plotNContentPlotly(fdl,
                         clusterNames = input$Ncluster,
                         usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    })
  }

  shinyApp(ui = ui, server = server)

}

