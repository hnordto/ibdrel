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

helpWidget = function(id) {
  tooltip(
    actionBttn(
      inputId = NS(id, "help"),
      label = NULL,
      style = "jelly",
      size = "sm",
      icon = icon("question"),
      color = "royal"
    ),
    title = "See help on this panel"
  )
}

# Load data XX FIX LATER XX

peds = ibdrel::ibdrel_unilineal$peds

annotation <- sapply(peds, annotatePedigree)
annotation.readable <- sapply(annotation, classTranslator, "eqclass.detailed")
names(annotation) <- annotation.readable
names(peds) <- annotation

segmentData = ibdrel::ibdrel_unilineal$segments
names(segmentData) <- annotation
metadata = ibdrel::pedsMetadata(peds)

supported.features <- c("count", "total", "median", "longest", "shortest")
input.types <- c("Segments" = "segments", "Summaries (not yet supported!)" = "summaries")

featureSelection.default = list("eqclass.detailed" = c("count", "length"),
                                "eqclass" = c("count", "length"),
                                "kappa" = c("count", "total"),
                                "kinship" = c("count", "total"),
                                "degree" = c("count", "total"))

inputBoxUI <- function(id) {

  tags$script(HTML("
  Shiny.addCustomMessageHandler('clickFile', function(message) {
    document.querySelector('[id$=\"segfile\"]').click();
  });
  "))

  bs4Card(
    title = div(class = "box-title-flex",
                div(class = "leftcolumn", "IBD Data"),
                div(class = "rightcolumn", helpWidget(id))
    ),
    status = "olive",
    width = NULL,
    collapse = TRUE,
    collapsed = FALSE,
    uiOutput(NS(id, "summaryField")),
    uiOutput(NS(id, "segmentsField")),

    # Hidden upload
    div(
      style = "display: none;",
      fileInput("segfile", NULL)
    ),

    footer = div(
      actionBttn(NS(id,"edit"), icon("edit"), style = "jelly", color = "primary", size = "m",
                 title = "Edit input"),
      actionBttn(NS(id,"load"), icon("folder-open"), style = "jelly", color = "primary", size = "m",
                 title = "Upload"),
      actionBttn(NS(id,"example"), icon("dice"), style = "jelly", color = "warning", size = "m",
                 title = "Simulate example")
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
      input_type = data[["inputType"]]

      resolutions <- c("eqclass.detailed", "eqclass", "kappa", "kinship", "degree")
      data[["resolutions"]] <- resolutions

      data[["segments"]] = sapply(resolutions, function(x) {
        aggregateSegments(segmentData, metadata, x)
      }, simplify = FALSE)

      data[["metadata"]] = metadata
      data[["peds"]] = peds

     # data[["features"]] <- lapply(data[["segments"]], function(x) {
    #    lapply(x, prepareFeatures, featureSel = featureSelection, cutoff = cutoff)
    #  })

      if (input_type == "segments") {
        features <- lapply(names(data[["segments"]]), function(nm) {
          x <- data[["segments"]][[nm]]
          lapply(x, prepareFeatures, featureSel = featureSelection.default[[nm]], cutoff = cutoff)
        })
        names(features) = names(data[["segments"]])
        data[["features"]] <- features



        #  data[["pdfs"]] <- lapply(data[["segments"]], function(x) {
        #    lapply(x, preparePdfs, featureSel = featureSelection, cutoff = cutoff)
        #  })

        pdfs <- lapply(names(data[["segments"]]), function(nm) {
          x <- data[["segments"]][[nm]]
          lapply(x, preparePdfs, featureSel = featureSelection.default[[nm]], cutoff = cutoff)
        })
        names(pdfs) <- names(data[["segments"]])
        data[["pdfs"]]  <- pdfs
      } else if (input_type == "summaries") {
        features <- lapply(data[["segments"]], function(x) {
          lapply(x, prepareFeatures, featureSel = featureSelection, cutoff = cutoff)
        })

        pdfs <- lapply(data[["segments"]], function(x) {
          lapply(x, preparePdfs, featureSel = featureSelection, cutoff = cutoff)
        })

      }

      # Features, only count+total for mdists
      ct <- lapply(data[["segments"]], function(x) {
        lapply(x, prepareFeatures, featureSel = c("count", "total"), cutoff = cutoff)
      })
      data[["ct"]] <- ct

    })

    # Change resolution
    observe({
      req(data[["selectedTab"]])
      data[["metadata"]]$class = data[["metadata"]][[data[["selectedTab"]]]]
      data[["classes"]] = names(data[["segments"]][[data[["selectedTab"]]]])
    })


    # Likelihood computation
    observe({
      req(!is.null(segment_input()))
      req(data[["cutoff"]])
      req(data[["outlierThreshold"]])
      req(data[["featureSelection"]])

      is_features <- data[["isFeatures"]]
      cutoff = data[["cutoff"]]
      input_type = data[["inputType"]]
      outlierThreshold = data[["outlierThreshold"]]
      selected_features = data[["featureSelection"]]

      if(input_type == "segments") {
        obs = as.numeric(segment_input() |>
                           strsplit("\n") |>
                           unlist())
        obs = obs[obs >= cutoff]
        data[["obs"]] = obs

        data[["likelihoods"]] = lapply(data[["pdfs"]], function(x) {
          classify(data[["obs"]], x, obsType = input_type, cutoff = cutoff, sort = FALSE) #!!!
        })

        data[["mdists"]] = lapply(data[["ct"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = c("count", "total"),
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "dist")
        })

        data[["mdistthreshold"]] = lapply(data[["ct"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = c("count", "total"),
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

        data[["mdists"]] = lapply(data[["ct"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = c("count", "total"),
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "dist")
        })

        data[["mdistthreshold"]] = lapply(data[["ct"]], function(x) {
          l <- distance(data[["obs"]], x, cutoff = cutoff, featureSel = c("count", "total"),
                        threshold = outlierThreshold, isFeatures = is_features)
          sapply(l, `[[`, "threshold")
        })
      }


    })

    observeEvent(input$edit, {
      input_type = data[["inputType"]]

      if (input_type == "segments") {
        showModal(modalDialog(
          uiOutput(NS(id, "segments_ui")),
          footer = tagList(
            actionButton(NS(id, "save_input_segments"), "Save")
          )
        ))
      } else if (input_type == "summaries") {
        showModal(modalDialog(
          uiOutput(NS(id, "summaries_ui"))
        ))
      }
    })


    observeEvent(input$load, {
      input_type = data[["inputType"]]

      if (input_type == "segments") {

        showModal(modalDialog(
          fileInput(NS(id, "segfile"), NULL,
                    width = "100%"),
          easyClose = TRUE
        ))
      }

    })

    observeEvent(input$segfile, {
      x = scan(input$segfile$datapath,
               what = numeric(),
               sep = ",")
      segment_input(paste(x, collapse = "\n"))
    })

    observeEvent(input$example, {
      showModal(modalDialog(
        selectInput(
        inputId = NS(id,"simrel"),
        width = "100%",
        label = "Relationship",
        choices = annotation
      ),
      footer = tagList(
        actionButton(
          inputId = NS(id,"simulateButton"),
          label = "Simulate",
          size = "sm",
          width = "100%",
          style = "simple",
          icon = icon("check")
          )
        )
      ))
    })
    observeEvent(input$simulateButton, {
      ped = data[["peds"]][[input$simrel]]
      simids = ibdrel::identifyLeaves(ped)

      segs = ibdsim2::ibdsim(ped, ids = simids, N = 1,
                             verbose = FALSE) |>
        ibdsim2::findPattern(pattern = list(carriers = simids)) |>
        ibdsim2::segmentStats(returnAll = TRUE)

      lens = segs$allSegs |> round(digits = 4)
      segment_input(paste(lens, collapse ="\n"))

      removeModal()
    })

    # Segment input UI

    output$summaryField <- renderUI({

      if(is.null(data[["obs"]])) return(tags$p("No data loaded."))

      segments.write = if (length(data[["obs"]]) == 1) "segment" else "segments"

      if(data[["inputType"]] == "segments") {
        tags$p(length(data[["obs"]]), segments.write, "above threshold", data[["cutoff"]], " cM with",
                             sum(data[["obs"]]), "cM total length.")
      } else if(data[["inputType"]] == "summaries") {
        count = if ("count" %in% names(data[["obs"]])) data[["obs"]][["count"]] else NULL
        total = if("total" %in% names(data[["obs"]])) data[["obs"]][["total"]] else NULL

        if(all(c("count", "total")  %in% names(data[["obs"]]))) {
          tags$p(count, "segments with", total, "cM total length.")
        } else if("count" %in% names(data[["obs"]])) {
          tags$p(count, "segments")
        } else if("total" %in% names(data[["obs"]])) {
          tags$p(total, "cM total length.")
        } else {
          tags$p("")
        }
      }
    })

    output$segmentsField <- renderUI({
      req(data[["obs"]])
      if(data[["inputType"]] != "segments") return(tags$p(""))

      tags$p("Individual segments:", paste(round(data[["obs"]],2), collapse = ", "))

    })


    # Remember segment input
    segment_input <- reactiveVal(NULL)
    observeEvent(input$save_input_segments, {
      segment_input(input$segmentInput)
      removeModal()
    }, ignoreInit = TRUE)

    output$segments_ui <- renderUI({
      req(data[["inputType"]])

      if (data[["inputType"]] != "segments") return(NULL)

      textAreaInput(
        inputId = NS(id,"segmentInput"),
        label = "Segment lengths (cM)",
        rows = 8,
        width = "100%",
        placeholder = "One segment per line, in cM",
        value = segment_input()
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

    # Help
    observeEvent(input$help, {
      shinyalert(
        className = "helpbox",
        html = TRUE,
        text = read_file("help/inputBox.html"),
        showConfirmButton = FALSE,
        closeOnClickOutside = TRUE,
        size = "m"
      )
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
    label = "Number of segments (count)",
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
