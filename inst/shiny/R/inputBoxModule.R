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

})


segmentDataRel = readRDS(system.file("data", "segments_unilineal_rel.rds", package = "ibdrel"))


# Tmp aggregate function: Implement in ibdrel later

aggregate.segments <- function(segments, metadata, metadata.agg.column) {
  aggregatedSegments = list()

  for (i in 1:length(segments)) {
    rel.id = names(segments)[i]
    agg.id = metadata |>
      dplyr::filter(rel == rel.id) |>
      slice(1) |>
      dplyr::select(metadata.agg.column) |>
      as.character()

    segments.rel <- segments[[i]]

    for (segment in segments.rel) {
      if (agg.id %in% names(aggregatedSegments)) {
        aggregatedSegments[[agg.id]][[length(aggregatedSegments[[agg.id]])+1]] <- segment
      } else {
        aggregatedSegments[[agg.id]][[1]] <- segment
      }
    }
  }

  return (aggregatedSegments)

}





#segmentDataKappa = readRDS(system.file("data", "segments_unilineal_kappa.rds", package = "ibdrel"))
#segmentDataKinship = readRDS(system.file("data", "segments_unilineal_kinship.rds", package = "ibdrel"))
segmentDataDeg = readRDS(system.file("data", "segments_unilineal_degree.rds", package = "ibdrel"))

pedsDataRel = readRDS(system.file("data", "peds_unilineal.rds", package = "ibdrel"))
metadata = ibdrel::pedsMetadata(pedsDataRel)

segmentData = readRDS(system.file("data", "segments_unilineal_rel.rds", package = "ibdrel"))
segmentData = aggregate.segments(segmentData,
                                 metadata,
                                 "eqclass.detailed")


inputBoxUI <- function(id) {

  box(
    width = NULL,
    collapsible = FALSE,
    title = "Input",

    # Model controls
    h5("Model"),
    pickerInput(inputId = NS(id, "classLevel"),
                label = "Resolution",
                choices = c("Sex-specific pedigree class" = "eqclass-detailed",
                            "Pedigree class" = "eqclass",
                            "Kappa" = "kappa",
                            "Kinship" = "kinship",
                            "Degree" = "degree"),
                multiple = FALSE),
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
    actionBttn(inputId = NS(id, "load"),
               label = "Load model"),

    # Input
    h5("Input"),
    tooltip(
      textAreaInput(inputId = NS(id, "segmentInput"),
                    label = "segment lengths",
                    rows = 10),
      title = "Segment input"
    )

  )

}

inputBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {
    observe({
      if (input$classLevel == "eqclass-detailed") {
        class.col = "eqclass.detailed"
        data[["segments"]] = segmentData # for now
      } else if (input$classLevel == "eqclass") {
        class.col = "eqclass"
        data[["segments"]] = segmentData
      } else if (input$classLevel == "degree") {
        class.col = "degree"
        data[["segments"]] = segmentDataDeg
      }

      metadata$class = metadata[[class.col]]
      data[["metadata"]] = metadata

      data[["classes"]] = names(data[["segments"]])

      updatePickerInput(session, "clasSelection",
                        choices = data[["classes"]],
                        selected = data[["classes"]])

    })

    observeEvent(input$load, {

      #data[["segments"]] <- data[["segments"]][names(data[["segments"]]) %in% input$clasSelection]


      data[["features"]] = lapply(data[["segments"]], prepareFeatures, cutoff = input$cutoff)
      data[["pdfs"]] = lapply(data[["segments"]], preparePdfs, cutoff = input$cutoff)

    })

    observeEvent(input$segmentInput, {
      obs = as.numeric(input$segmentInput |>
                         strsplit("\n") |>
                         unlist())
      obs = obs[obs > input$cutoff]
      data[["obs"]] = obs
    })

    observeEvent(input$segmentInput, {
      req(data[["obs"]], data[["pdfs"]])
      likelihoods = classify(data[["obs"]], data[["pdfs"]], input$cutoff, sort = TRUE)
      data[["likelihoods"]][["raw"]] = likelihoods
      data[["likelihoods"]][["norm"]] = ibdrel::normalizePosteriors(likelihoods)

      #data[["mdists"]] = distance(data[["obs"]], data[["features"]], orderLst = data[["likelihoods"]], top_n = input$n_outlier)
    })
  })

}
