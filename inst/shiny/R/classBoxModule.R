checkInput = function(data) {
  if (is.null(data)) {
    validate(need(FALSE, "Waiting for input."))
  }
  data
}

classBoxUI = function(id) {

  tabBox(
    id = NS(id, "classTabs"),
    width = NULL,
    selected = "Summary",

    title = div(
      pickerInput(inputId = NS(id, "selectClass"), label = "Class:", choices = NULL)
    ),

    tabPanel(
      title = "Summary",
      tags$h5("Summary"),
      uiOutput(NS(id, "classInfo")),
      tags$h5("IBD segment distribution"),
      plotOutput(NS(id, "jointDistPlot")),
      tags$h5("Outlier diagnostics"),
      uiOutput(NS(id, "outlierInfo"))
    ),

    tabPanel(
      title = "Pedigrees",
      div(
        class = "table-box",
        gt_output(NS(id, "relationshipTable"))
      ),
      selectInput(NS(id, "selectPedigree"), "Pedigree:", choices = NULL)
    )

  )
}

classBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {


    #observe({
    #  updateSelectInput(session, "selectPedigree",
    #                    choices = subset(data[["metadata"]], class == input$selectClass)[["rel"]])
    #})

    observe({
      updatePickerInput(session, "selectClass", choices = data[["classes"]])
      #data[["selectedClass"]] = input$selectClass
    })

    output$classInfo <- renderUI({
      req(input$selectClass)

      likelihoods <- checkInput(data[["likelihoods"]])

      makeSummaryOutput(data[["metadata"]],
                        likelihoods,
                        input$selectClass,
                        resolution = data[["selectedTab"]])
    })


    output$outlierInfo <- renderUI({

      mdists <- checkInput(data[["mdists"]])

      makeOutlierOutput(mdists,
                        input$selectClass,
                        data[["selectedTab"]],
                        data[["mdistthreshold"]])
    })

    output$relationshipTable = render_gt({
      relationshipTable(data[["metadata"]], input$selectClass)
    })

    output$jointDistPlot = renderPlot({

      resolution = data[["selectedTab"]]

      if (!is.null(data[["obs"]])) {
        jointDistPlot(train.features = data[["features"]][[resolution]][[input$selectClass]],
                      class = input$selectClass,
                      obs.features = ibdrel::obsToFeatures(data[["obs"]]))
      } else {
        jointDistPlot(train.features = data[["features"]][[resolution]][[input$selectClass]],
                      class = input$selectClass,
                      obs.features = NULL)
      }
    })


  })
}

makeSummaryOutput <- function(metadata, likelihoods, selected.class, resolution) {

  met = subset(metadata, class == selected.class)

  tagList(
    tags$p(
      tags$b("Selected class:"), tags$i(ibdrel::classTranslator(selected.class, resolution)), "\n",
      tags$ul(
        tags$li(tags$b("Log-likelihood:"), likelihoods[["raw"]][[resolution]][[selected.class]]),
        tags$li(tags$b("Normalized likelihood:"), likelihoods[["norm"]][[resolution]][[selected.class]]),
        tags$li(tags$b("Rank:"), which(names(likelihoods[["raw"]][[resolution]]) == selected.class)),
        tags$li(tags$b("Number of pedigrees in class:"), length(met[["code"]]))
      )
    )
  )
}

makeOutlierOutput <- function(mdists, selected.class, resolution, threshold) {

  mdist.idx = which(names(mdists[[resolution]]) == selected.class)
  mdist = mdists[[resolution]][mdist.idx]

  tagList(
    tags$p(
      tags$b("Is observation an outlier?"), ifelse(mdist <= threshold, "No", "Yes"), "\n",
      tags$ul(
        tags$li(tags$b("Mahalanobis distance:"), mdist),
        tags$li(tags$b("Threshold:"), threshold)
      )
    )
  )
}

relationshipTable <- function(metadata, selected.class) {
  df <- subset(metadata, class == selected.class)
  df <- df[,c("code", "kappa0", "kappa1", "kappa2", "kinship")]
  df$kappa0 <- as.character(MASS::fractions(df$kappa0))
  df$kappa1 <- as.character(MASS::fractions(df$kappa1))
  df$kappa2 <- as.character(MASS::fractions(df$kappa2))
  df$kinship <- as.character(MASS::fractions(df$kinship))
  df <- gt(df) |>
    gt_theme_538() #|>
    #cols_label(
    #  kappa0 = html("&kappa;<sub>0</sub>"),
    #  kappa1 = html("&kappa;<sub>1</sub>"),
    #  kappa2 = html("&kappa;<sub>2</sub>"),
     # kinship = html("&phi;")
    #)
  return (df)
}


jointDistPlot <- function(train.features,
                          class,
                          obs.features = NULL) {
  count = train.features[["countPdf"]]
  total = train.features[["totalPdf"]]

  data = data.frame(count, total)

  p <- ggplot2::ggplot() +
    geom_jitter(data = data, mapping = aes(x = count, y = total), alpha = .5) +
    #stat_ellipse(data = data, mapping = aes(x = count, y = total, size = 1.1), alpha = .75) +
    labs(x = "Number of segments",
         y = "Total IBD segment length (cM)") +
    theme_bw()

  if (!is.null(obs.features)) {
    count.obs = obs.features$obs[["countPdf"]]
    total.obs = obs.features$obs[["totalPdf"]]

    data.obs = data.frame(count = count.obs, total = total.obs)

    p <- p +
      geom_point(data = data.obs, mapping = aes(x = count, y = total), alpha = .5,
                  shape = 3, stroke = 2, size = 4, colour = "red2")
  }

  return (p)

}
