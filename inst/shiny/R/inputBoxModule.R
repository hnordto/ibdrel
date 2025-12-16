suppressPackageStartupMessages({
  library(shiny)
  library(ibdrel)
  library(ibdsim2)
  library(bs4Dash)
  library(ggplot2)
  library(stringr)
  library(gt)
  library(DT)
  library(scales)
  library(dplyr)
  library(gtExtras)
  library(shinyWidgets)
  library(pedtools)
  library(shinyalert)
  library(readr)
  library(shinyjs)

})

# Load data XX FIX LATER XX

peds = ibdrel::ibdrel_unilineal$peds

annotation <- sapply(peds, annotatePedigree)
names(peds) <- annotation

segmentData = ibdrel::ibdrel_unilineal$segments
names(segmentData) <- annotation
metadata = ibdrel::pedsMetadata(peds)

supported.features <- c("count", "total", "median", "longest", "shortest")

inputBoxUI <- function(id) {

  box(
    width = NULL,
    collapsible = FALSE,
    title = "Input",

    # Model controls
    h5("Model"),
    pickerInput(inputId = NS(id, "featureSelection"),
                label = "Features",
                choices = supported.features,
                selected = c("count", "total"),
                options = pickerOptions(
                  actionsBox = TRUE,
                  size = 10,
                  selectedTextFormat = "count > 3"
                ),
                multiple = TRUE),
    numericInput(inputId = NS(id, "cutoff"),
                 label = "Cutoff",
                 value = 7, min = 0, step = 1),

    # Input
    h5("Input"),
    tabsetPanel(
      id = NS(id, "inputMode"),
      type = "pills",

      tabPanel(
        title = "Segments",
        value = "segments"
      ),

      tabPanel(
        title = "Simulate",
        value = "simulate",
        tooltip(
          selectInput(inputId = NS(id, "simrel"),
                      label = "Relationship",
                      choices = annotation),
          title = "Simulate"),
        actionBttn(inputId = NS(id, "simulateButton"),
                   "Simulate")
      ),

      tabPanel(
        title = "Upload",
        value = "upload",
        tooltip(
          fileInput(inputId = NS(id, "segfile"),
                    label = "Segments file",
                    buttonLabel = icon("folder-open"),
                    accept = c(".txt")),
          title = "Upload"
        )

      )

    ),

    tooltip(
      textAreaInput(inputId = NS(id, "segmentInput"),
                    label = "Segment lengths (cM)",
                    rows = 10),
      title = "Segment input"
    ),

    # Controls
    h5("Results"),
    materialSwitch(inputId = NS(id, "normalizeButton"), value = TRUE,
                   label = "Normalize likelihoods"),
    materialSwitch(inputId = NS(id, "filterButton"), value = FALSE,
                   label = "Filter outliers"),

    h5("Outlier detection"),
    numericInput(
      inputId = NS(id, "outlierThreshold"),
      label = HTML("&#967;<sup>2</sup> quantile threshold"),
      min = 0,
      max = 1,
      value = 0.95,
      step = 0.01
    )

  )

}

inputBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {
    observe({

      resolutions <- c("eqclass.detailed", "eqclass", "kappa", "kinship", "degree")
      data[["resolutions"]] <- resolutions

      data[["segments"]] = sapply(resolutions, function(x) {
        aggregateSegments(segmentData, metadata, x)
      }, simplify = FALSE)

      data[["metadata"]] = metadata
      data[["peds"]] = peds

      data[["features"]] <- input$featureSelection

      data[["features"]] <- lapply(data[["segments"]], function(x) {
        lapply(x, prepareFeatures, featureSel = input$featureSelection, cutoff = input$cutoff)
      })

      data[["pdfs"]] <- lapply(data[["segments"]], function(x) {
        lapply(x, preparePdfs, featureSel = input$featureSelection, cutoff = input$cutoff)
      })


    })

    observe({
      req(data[["selectedTab"]])
      data[["metadata"]]$class = data[["metadata"]][[data[["selectedTab"]]]]
      data[["classes"]] = names(data[["segments"]][[data[["selectedTab"]]]])
    })



    observeEvent(input$segmentInput, {

      obs = as.numeric(input$segmentInput |>
                         strsplit("\n") |>
                         unlist())
      obs = obs[obs >= input$cutoff]
      data[["obs"]] = obs

      data[["mdists"]] = lapply(data[["features"]], function(x) {
        distance(data[["obs"]], x, cutoff = input$cutoff, featureSel = input$featureSelection)
      })
    })

    observe({
      data[["likelihoods"]] = lapply(data[["pdfs"]], function(x) {
        classify(data[["obs"]], x, cutoff = input$cutoff, sort = FALSE) #!!!
      })
    })


    # Controls
    observe({
      data[["normalize"]] = input$normalizeButton
      data[["chisq"]] = input$outlierThreshold
      data[["mdistthreshold"]] = qchisq(p = input$outlierThreshold,
                                       df = length(data[["features"]][[1]][[1]])) # Arbitrary feature vector
      data[["filter"]] = ifelse(input$filterButton, TRUE, FALSE)
      data[["cutoff"]] = input$cutoff
    })

    # Input mode

    observeEvent(input$inputMode, {
      if (input$inputMode %in% c("simulate", "upload")) {
        shinyjs::runjs(sprintf(
          "$('#%s').prop('readonly', true);",
          NS(id, "segmentInput")
        ))
      } else {
        shinyjs::runjs(sprintf(
          "$('#%s').prop('readonly', false);",
          NS(id, "segmentInput")
        ))
      }
    }, ignoreInit = TRUE)


    # Simulate example

    observeEvent(input$simulateButton, {

      ped = data[["peds"]][[input$simrel]]
      simids = ibdrel::identifyLeaves(ped)

      segs = ibdsim2::ibdsim(ped, ids = simids, N = 1,
                             verbose = FALSE) |>
        ibdsim2::findPattern(pattern = list(carriers = simids)) |>
        ibdsim2::segmentStats(returnAll = TRUE)

      lens = segs$allSegs |> round(digits = 4)

      updateTextAreaInput(session, "segmentInput", value = paste(lens, collapse = "\n"))

    })

    # Load from file

    observeEvent(input$segfile, {
      x = scan(input$segfile$datapath,
               what = numeric(),
               sep = ",")
      updateTextAreaInput(session, "segmentInput", value = paste(x, collapse = "\n"))
    })
  })

}


