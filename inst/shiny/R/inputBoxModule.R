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

# Data
segmentData = readRDS(system.file("data", "segments_unilineal_rel.rds", package = "ibdrel"))
pedsDataRel = readRDS(system.file("data", "peds_unilineal.rds", package = "ibdrel"))
metadata = ibdrel::pedsMetadata(pedsDataRel)





#segmentDataKappa = readRDS(system.file("data", "segments_unilineal_kappa.rds", package = "ibdrel"))
#segmentDataKinship = readRDS(system.file("data", "segments_unilineal_kinship.rds", package = "ibdrel"))
#segmentDataDeg = readRDS(system.file("data", "segments_unilineal_degree.rds", package = "ibdrel"))



#segmentData = readRDS(system.file("data", "segments_unilineal_rel.rds", package = "ibdrel"))
#segmentData = aggregate.segments(segmentData,
#                                 metadata,
#                                 "eqclass.detailed")


inputBoxUI <- function(id) {

  box(
    width = NULL,
    collapsible = FALSE,
    title = "Input",

    # Model controls
    h5("Model"),
    pickerInput(inputId = NS(id, "clasSelection"),
                label = "Classes",
                choices = NULL,
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
    tooltip(
      textAreaInput(inputId = NS(id, "segmentInput"),
                    label = "Segment lengths (cM)",
                    rows = 10),
      title = "Segment input"
    ),

    # Controls
    h5("Results"),
    materialSwitch(inputId = NS(id, "normalizeButton"), value = TRUE,
                   label = "Normalize likelihoods?"),
    materialSwitch(inputId = NS(id, "filterButton"), value = FALSE,
                   label = "Filter outliers?"),

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

      #data[["classes"]] = lapply(data[["segments"]], function(x) {
      #  names(x)
      #})

      data[["features"]] <- lapply(data[["segments"]], function(x) {
        lapply(x, prepareFeatures, cutoff = input$cutoff)
      })

      data[["pdfs"]] <- lapply(data[["segments"]], function(x) {
        lapply(x, preparePdfs, cutoff = input$cutoff)
      })


    })

    observe({
      req(data[["selectedTab"]])
      data[["metadata"]]$class = data[["metadata"]][[data[["selectedTab"]]]]
      data[["classes"]] = names(data[["segments"]][[data[["selectedTab"]]]])

      updatePickerInput(session, "clasSelection",
                        choices = data[["classes"]],
                        selected = data[["classes"]])
    })



    observeEvent(input$segmentInput, {

      obs = as.numeric(input$segmentInput |>
                         strsplit("\n") |>
                         unlist())
      obs = obs[obs > input$cutoff]
      data[["obs"]] = obs

      data[["mdists"]] = lapply(data[["features"]], function(x) {
        distance(data[["obs"]], x)
      })
    })

    observeEvent(input$filterButton, {

      data[["likelihoods"]][["raw"]] = lapply(data[["pdfs"]], function(x) {
        classify(data[["obs"]], x, cutoff = input$cutoff, sort = TRUE)
      })
      data[["likelihoods"]][["norm"]] = lapply(data[["likelihoods"]][["raw"]],
                                               ibdrel::normalizeLikelihoods)

      if (input$filterButton) {
        data[["likelihoods"]][["raw"]] = Map(
          function(x, mdists) {
            mdists.ordered <- mdists[names(x)]
            x[mdists.ordered <= data[["mdistthreshold"]]]
          },
          data[["likelihoods"]][["raw"]], data[["mdists"]]
        )

        data[["likelihoods"]][["norm"]] = lapply(data[["likelihoods"]][["raw"]],
                                                 ibdrel::normalizeLikelihoods)

      }
    })

    # Controls
    observe({
      data[["normalize"]] = input$normalizeButton
      data[["mdistthreshold"]] = qchisq(p = input$outlierThreshold,
                                        df = length(data[["features"]][[1]][[1]])) # Arbitrary feature vector
    })


  })

}
