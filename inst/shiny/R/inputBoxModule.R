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
input.types <- c("Segments" = "segments", "Summaries" = "summaries")

inputBoxUI <- function(id) {

  bs4Card(
    title = "IBD Data",
    width = NULL,
    collapse = TRUE,
    collapsed = FALSE,
    verbatimTextOutput(NS(id, "summaryField")),
    verbatimTextOutput(NS(id, "segmentsField")),
    footer = div(
      actionBttn(NS(id,"edit"), icon("edit"), style = "jelly", color = "primary", size = "m"),
      actionBttn(NS(id,"load"), icon("folder-open"), style = "jelly", color = "primary", size = "m"),
      actionBttn(NS(id,"example"), icon("dice"), style = "jelly", color = "warning", size = "m")
    )
  )

  #box(
  #  width = NULL,
  #  collapsible = FALSE,
  #  title = NULL,
  #  class = "input-box",


    #conditionalPanel(
    #  condition =  sprintf("input['%s'] == 'segments'", NS(id,"inputType")),
    #  h5("Segments"),
    #  tooltip(
    #    div(
    #      class = "file-upload",
    #      fileInput(
    #        inputId = NS(id,"segfile"),
    #        label = "Upload segments file (.txt)",
    #        buttonLabel = icon("folder-open"),
    #        accept = ".txt"
    #      )
    #    ),
    #    title = "Upload file with segment lengths"
    #  ),

    #  tags$label("Simulate example"),
    #  shinyWidgets::dropdown(
    #    inputId = NS(id, "simulateDropdown"),
    #    label = "Simulate",
    #    icon = icon("dice"),
    #    circle = FALSE,
    #    status = "default",
    #    width = "100%",
    #    tooltip(
    #      selectInput(
    #        inputId = NS(id,"simrel"),
    #        label = "Relationship",
    #        choices = annotation
    #      ),
    #      title = "Relationship used for simulation"
    #    ),
    #    actionButton(
    #      inputId = NS(id,"simulateButton"),
    #      label = "Simulate",
    #      size = "sm",
    #      width = "100%",
    #      style = "simple",
    #      icon = icon("check")
    #   )
    #  )
    #),

    #tagList(
    #  uiOutput(NS(id, "segments_ui"))
   # ),

  #  tagList(
   #   uiOutput(NS(id, "summaries_ui"))
  #  ),

    #conditionalPanel(
    #  condition =  sprintf("input['%s'] == 'summaries'", NS(id,"inputType")),
    #  uiOutput(NS(id, "summaries_ui")),
    #),



   # hr()
    # use dropdown() instead of dropdownButton() because the latter don't work with pickerInput()
    #shinyWidgets::dropdown(
    #  inputId = NS(id, "settingsDropdown"),
    #  label = "Settings",
    #  icon = icon("gear"),
    #  circle = FALSE,
    #  size = "sm",
    #  status = "default",
    #  width = "100%",
    #  pickerInput(inputId = NS(id, "featureSelection"),
    #              label = "Features",
    #              choices = supported.features,
    #              selected = c("count", "total"),
    #              options = pickerOptions(
    #                actionsBox = TRUE,
    #                size = 10,
    #                selectedTextFormat = "count > 3"
    #              ),
    #              multiple = TRUE),
    #  materialSwitch(inputId = NS(id, "normalizeButton"), value = TRUE,
    #                 label = "Normalize likelihoods"),
    #  materialSwitch(inputId = NS(id, "filterButton"), value = FALSE,
    #                 label = "Filter outliers"),

    #  h5("Outlier detection"),
    #  numericInput(
    #    inputId = NS(id, "outlierThreshold"),
    #    label = HTML("Outlier threshold"),
    #    min = 0,
    #    max = 1,
    #    value = 0.95,
    #    step = 0.01
    #  )
    #)

  #)

}

inputBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    ns <- session$id

    # Initialize model
    observe({

      req(data[["cutoff"]])

      cutoff = data[["cutoff"]]
      featureSelection = data[["featureSelection"]]

      resolutions <- c("eqclass.detailed", "eqclass", "kappa", "kinship", "degree")
      data[["resolutions"]] <- resolutions

      data[["segments"]] = sapply(resolutions, function(x) {
        aggregateSegments(segmentData, metadata, x)
      }, simplify = FALSE)

      data[["metadata"]] = metadata
      data[["peds"]] = peds

      data[["features"]] <- lapply(data[["segments"]], function(x) {
        lapply(x, prepareFeatures, featureSel = featureSelection, cutoff = cutoff)
      })

      data[["pdfs"]] <- lapply(data[["segments"]], function(x) {
        lapply(x, preparePdfs, featureSel = featureSelection, cutoff = cutoff)
      })


    })

    # Change resolution
    observe({
      req(data[["selectedTab"]])
      data[["metadata"]]$class = data[["metadata"]][[data[["selectedTab"]]]]
      data[["classes"]] = names(data[["segments"]][[data[["selectedTab"]]]])
    })


    # Likelihood computation
    observe({
      req(!is.null(input$segmentInput))
      req(data[["cutoff"]])
      req(data[["outlierThreshold"]])
      req(data[["featureSelection"]])

      is_features <- data[["isFeatures"]]
      cutoff = data[["cutoff"]]
      input_type = data[["inputType"]]
      outlierThreshold = data[["outlierThreshold"]]
      selected_features = data[["featureSelection"]]

      if(input_type == "segments") {
        obs = as.numeric(input$segmentInput |>
                           strsplit("\n") |>
                           unlist())
        obs = obs[obs >= cutoff]
        data[["obs"]] = obs

        data[["likelihoods"]] = lapply(data[["pdfs"]], function(x) {
          classify(data[["obs"]], x, obsType = input_type, cutoff = cutoff, sort = FALSE) #!!!
        })

        data[["mdists"]] = lapply(data[["features"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = selected_features,
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "dist")
        })

        data[["mdistthreshold"]] = lapply(data[["features"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = selected_features,
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "threshold")
        })
      } else if(input_type == "summaries") {
        values <- lapply(selected_features, function(x) {
          input[[x]]
        })

        req(all(!sapply(values, is.null)))
        req(all(!is.na(values)))

        summaries <- setNames(values, selected_features)

        data[["obs"]] <- summaries

        data[["likelihoods"]] = lapply(data[["pdfs"]], function(x) {
          classify(data[["obs"]], x, obsType = input_type, cutoff = cutoff, sort = FALSE) #!!!
        })

        data[["mdists"]] = lapply(data[["features"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = selected_features,
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "dist")
        })

        data[["mdistthreshold"]] = lapply(data[["features"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = selected_features,
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "threshold")
        })
      }


    })

    observeEvent(input$edit, {
      input_type = data[["inputType"]]
      showModal(modalDialog(
        if(input_type == "segments") {
          tagList(
            uiOutput(NS(id, "segments_ui"))
          )
        } else if(input_type == "summaries") {
          tagList(
            uiOutput(NS(id, "summaries_ui"))
          )
        }
      ))
    })

    observeEvent(input$example, {
      showModal(modalDialog(
              selectInput(
                inputId = NS(id,"simrel"),
                label = "Relationship",
                choices = annotation
              ),
              actionButton(
                inputId = NS(id,"simulateButton"),
                label = "Simulate",
                size = "sm",
                width = "100%",
                style = "simple",
                icon = icon("check")
             )
      ))
    })



    # Simulate example


    observeEvent(input$simulateButton, {

      ped = data[["peds"]][[input$simrel]]
      simids = ibdrel::identifyLeaves(ped)

      segs = ibdsim2::ibdsim(ped, ids = simids, N = 1,
                             verbose = FALSE) |>
        ibdsim2::findPattern(pattern = list(carriers = simids)) |>
        ibdsim2::segmentStats(returnAll = TRUE)

      lens = segs$allSegs |> round(digits = 4)

      updateTextAreaInput(session, NS(id,"segmentInput"), value = paste(lens, collapse = "\n"))

    })

    # Load from file

    observeEvent(input$segfile, {
      x = scan(input$segfile$datapath,
               what = numeric(),
               sep = ",")
      updateTextAreaInput(session, "segmentInput", value = paste(x, collapse = "\n"))
    })

    # Segment input UI

    output$summaryField <- renderPrint({

      if(is.null(data[["obs"]])) return(cat("No data loaded."))

      if(data[["inputType"]] == "segments") {
        cat(length(data[["obs"]]), "segments with",
                             sum(data[["obs"]]), "cM total length.")
      } else if(data[["inputType"]] == "summaries") {
        count = if ("count" %in% names(data[["obs"]])) data[["obs"]][["count"]] else NULL
        total = if("total" %in% names(data[["obs"]])) data[["obs"]][["total"]] else NULL

        if(all(c("count", "total")  %in% names(data[["obs"]]))) {
          cat(count, "segments with", total, "cM total length.")
        } else if("count" %in% names(data[["obs"]])) {
          cat(count, "segments")
        } else if("total" %in% names(data[["obs"]])) {
          cat(total, "cM total length.")
        } else {
          cat(NULL)
        }
      }
    })

    output$segmentsField <- renderPrint({
      req(data[["obs"]])
      if(data[["inputType"]] != "segments") return(NULL)

      cat("Individual segments:", paste(data[["obs"]], collapse = ","))

    })



    output$segments_ui <- renderUI({
      req(data[["inputType"]])

      if (data[["inputType"]] != "segments") return(NULL)

      textAreaInput(
        inputId = NS(id,"segmentInput"),
        label = "Segment lengths (cM)",
        rows = 8,
        placeholder = "One segment per line, in cM"
      )
    })


    # Summaries UI

    output$summaries_ui <- renderUI({
      req(data[["featureSelection"]])
      req(data[["inputType"]])

      if (data[["inputType"]] != "summaries") return(NULL)

      selected_features <- data[["featureSelection"]]

      ui_list <- lapply(selected_features, function(feat) {
        spec <- feature_specs[[feat]]

        numericInput(
          inputId = NS(id,paste0(feat)),
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

summaries_obs <- list(
  count = integer(),
  total = numeric(),
  median = numeric(),
  longest = numeric(),
  shortest = numeric()
)

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
