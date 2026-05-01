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
    title = NULL,

    h5("IBD input"),
    radioButtons(
      inputId = NS(id, "inputType"),
      label = NULL,
      inline = TRUE,
      choices = c(
        "Segments" = "segments",
        "Summaries" = "summaries"
      )
    ),

    conditionalPanel(
      condition =  sprintf("input['%s'] == 'segments'", NS(id,"inputType")),
      h5("Segments"),
      tooltip(
        textAreaInput(
          inputId = NS(id,"segmentInput"),
          label = "Segment lengths (cM)",
          rows = 8,
          placeholder = "One segment per line, in cM"
        ),
        title = "Paste or type segment lengths"
      ),
      tooltip(
        div(
          class = "file-upload",
          fileInput(
            inputId = NS(id,"segfile"),
            label = "Upload segments file (.txt)",
            buttonLabel = icon("folder-open"),
            accept = ".txt"
          )
        ),
        title = "Upload file with segment lengths"
      ),

      tags$label("Simulate example"),
      shinyWidgets::dropdownButton(
        inputId = NS(id, "simulateDropdown"),
        label = "Simulate",
        icon = icon("dice"),
        circle = FALSE,
        status = "default",
        width = "100%",
        tooltip(
          selectInput(
            inputId = NS(id,"simrel"),
            label = "Relationship",
            choices = annotation
          ),
          title = "Relationship used for simulation"
        ),
        actionButton(
          inputId = NS(id,"simulateButton"),
          label = "Simulate",
          size = "sm",
          style = "simple",
          width = "100%",
          icon = icon("check")
        )
      )
    ),

    conditionalPanel(
      condition =  sprintf("input['%s'] == 'summaries'", NS(id,"inputType")),
      uiOutput(NS(id, "summaries_ui"))

    ),

    hr(),
    h5("Match threshold"),
    numericInput(inputId = NS(id, "cutoff"),
                 label = NULL,
                 value = 7, min = 0, step = 1),
    hr(),
    # use dropdown() instead of dropdownButton() because the latter don't work with pickerInput()
    shinyWidgets::dropdown(
      inputId = NS(id, "settingsDropdown"),
      label = "Settings",
      icon = icon("gear"),
      circle = FALSE,
      status = "default",
      width = "100%",
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
      materialSwitch(inputId = NS(id, "normalizeButton"), value = TRUE,
                     label = "Normalize likelihoods"),
      materialSwitch(inputId = NS(id, "filterButton"), value = FALSE,
                     label = "Filter outliers"),

      h5("Outlier detection"),
      numericInput(
        inputId = NS(id, "outlierThreshold"),
        label = HTML("Outlier threshold"),
        min = 0,
        max = 1,
        value = 0.95,
        step = 0.01
      )
    )

  )

}

inputBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    ns <- session$id

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


      shinyjs::reset(NS(id,"segfile"))


      obs = as.numeric(input$segmentInput |>
                         strsplit("\n") |>
                         unlist())
      obs = obs[obs >= input$cutoff]
      data[["obs"]] = obs

      data[["mdists"]] = lapply(data[["features"]], function(x) {
        l <- distance(data[["obs"]], x, cutoff = input$cutoff, featureSel = input$featureSelection,
                      threshold = input$outlierThreshold)
        sapply(l, `[[`, "dist")
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
      data[["mdistthreshold"]] = lapply(data[["features"]], function(x) {
        l <- distance(data[["obs"]], x, cutoff = input$cutoff, featureSel = input$featureSelection,
                      threshold = input$outlierThreshold)
        sapply(l, `[[`, "threshold")
      })
      #data[["mdistthreshold"]] = qchisq(p = input$outlierThreshold,
      #                                 df = length(data[["features"]][[1]][[1]])) # Arbitrary feature vector
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

    # Summaries UI

    output$summaries_ui <- renderUI({
      req(input$featureSelection)
      req(input$inputType)

      if (input$inputType != "summaries") return(NULL)

      selected_features <- input$featureSelection

      ui_list <- lapply(selected_features, function(feat) {
        spec <- feature_specs[[feat]]

        numericInput(
          inputId = NS(id,paste0("summary_", feat)),
          label = spec$label,
          value = spec$value,
          min = spec$min,
          step = spec$step
        )
      })

      tagList(ui_list)

    })

  })
}


feature_specs <- list(
  count = list(
    label = "Numer of segments (count)",
    value = NA,
    min = 0,
    step = 1
  ),
  total = list(
    label = "Total segment length (cM)",
    value = NA,
    min = 0,
    step = 0.1
  ),
  median = list(
    label = "Median segment length (cM)",
    value = NA,
    min = 0,
    step = 0.1
  ),
  longest = list(
    label = "Longest segment length (cM)",
    value = NA,
    min = 0,
    step = 0.1
  ),
  shortest = list(
    label = "Shortest segment length (cM)",
    value = NA,
    min = 0,
    step = 0.1
  )
)
