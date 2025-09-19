checkLikelihood = function(likelihood) {
  if (all(likelihood == -Inf)) {
    validate(need(FALSE, "Waiting for input.\nNeed at least one segment with length >= cutoff."))
  }
}

checkObs = function(obs) {
  if(identical(obs, numeric(0))) {
    validate(need(FALSE, "Waiting for input.\nNeed at least one segment with length >= cutoff."))
  }
}

classBoxUI = function(id) {

  tabBox(
    id = NS(id, "classTabs"),
    width = NULL,
    selected = "Summary",
    collapsible = FALSE,

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
        class = "table-box-null",
        gt_output(NS(id, "relationshipTable"))
      ),
      pickerInput(NS(id, "selectPedigree"), "Pedigree:", choices = NULL),
      plotOutput(NS(id, "pedPlot"))
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

      res = data[["selectedTab"]]
      likelihoods = data[["likelihoods"]][[res]]

      # If not input, no ordering
      if (!is.null(names(likelihoods))) {
        classes <- names(sort(likelihoods, decreasing = TRUE))
      } else {
        classes <- data[["classes"]]
      }

      if(length(classes) == 0) return() # Escape if classes is not loaded yet (Shiny timing issue)

      classes.labels <- sapply(classes, ibdrel::classTranslator, res)
      updatePickerInput(session, "selectClass",
                        choices = stats::setNames(classes, classes.labels))

    })

    observeEvent(input$selectClass, {
      peds = subset(data[["metadata"]], class == input$selectClass)[["rel"]]
      updatePickerInput(session, "selectPedigree", choices = peds)
    })

    output$classInfo <- renderUI({
      req(input$selectClass)

      checkLikelihood(data[["likelihoods"]][[data[["selectedTab"]]]])

      makeSummaryOutput(data, data[["selectedTab"]], input$selectClass)
    })


    output$outlierInfo <- renderUI({
      req(input$selectClass)

      checkObs(data[["obs"]])

      makeOutlierOutput(data, data[["selectedTab"]], input$selectClass)
    })

    output$relationshipTable = render_gt({
      req(input$selectClass)
      relationshipTable(data, input$selectClass)
    })

    output$jointDistPlot = renderPlot({
      req(input$selectClass)

      resolution = data[["selectedTab"]]

      jointDistPlot(data, resolution, input$selectClass)
    })

    output$pedPlot <- renderPlot({
      req(input$selectPedigree)

      ped = data[["peds"]][[input$selectPedigree]]
      ped.leaves = ibdrel::identifyLeaves(ped)
      ped = ped |> pedtools::setSex(leaves, sex= 0)

      plot(ped, autoScale = TRUE, hatched = ped.leaves, fill = list(red = ped.leaves))

    })


  })
}

makeSummaryOutput <- function(data, resolution, selected.class) {

  metadata = data[["metadata"]]
  met = subset(metadata, class == selected.class)

  likelihoods = data[["likelihoods"]][[resolution]]
  likelihoods.sorted = sort(likelihoods, decreasing = TRUE)

  likelihood = likelihoods[[selected.class]]
  rank = which(names(likelihoods.sorted) == selected.class)
  npeds = length(met[["code"]])



  tagList(
    tags$p(
      tags$b("Selected class:"), tags$i(ibdrel::classTranslator(selected.class, resolution)), "\n",
      tags$ul(
        tags$li(tags$b("Log-likelihood:"), round(likelihood, 3)),
        #tags$li(tags$b("Normalized likelihood:"), likelihoods[["norm"]][[resolution]][[selected.class]]),
        tags$li(tags$b("Rank:"), rank),
        tags$li(tags$b("Number of pedigrees in class:"), npeds)
      )
    )
  )
}

makeOutlierOutput <- function(data, resolution, selected.class) {

  mdists = data[["mdists"]]
  threshold = data[["mdistthreshold"]]

  mdist.idx = which(names(mdists[[resolution]]) == selected.class)
  mdist = mdists[[resolution]][mdist.idx]

  tagList(
    tags$p(
      tags$b("Is observation an outlier?"), ifelse(mdist <= threshold, "No", "Yes"), "\n",
      tags$ul(
        tags$li(tags$b("Mahalanobis distance:"), round(mdist, 3)),
        tags$li(tags$b("Threshold:"), round(threshold, 3))
      )
    )
  )
}

relationshipTable <- function(data, selected.class) {

  metadata <- data[["metadata"]]


  df <- subset(metadata, class == selected.class)
  df <- df[,c("code", "kappa0", "kappa1", "kappa2", "kinship")]
  df$kappa0 <- as.character(MASS::fractions(df$kappa0))
  df$kappa1 <- as.character(MASS::fractions(df$kappa1))
  df$kappa2 <- as.character(MASS::fractions(df$kappa2))
  df$kinship <- as.character(MASS::fractions(df$kinship))
  df <- gt(df) |>
    gt_theme_538(quiet = TRUE) #|>
    #cols_label(
    #  kappa0 = html("&kappa;<sub>0</sub>"),
    #  kappa1 = html("&kappa;<sub>1</sub>"),
    #  kappa2 = html("&kappa;<sub>2</sub>"),
     # kinship = html("&phi;")
    #)
  return (df)
}


jointDistPlot <- function(data, resolution, selected.class) {

  segments = data[["features"]][[resolution]][[selected.class]]

  count = segments[["count"]]
  total = segments[["total"]]

  obs = data[["obs"]]

  if (!identical(obs, numeric(0))) {
    obs.features = ibdrel::obsToFeatures(obs, data[["cutoff"]], c("count", "total"))
  } else {
    obs.features = NULL
  }
  df = data.frame(count, total)

  p <- ggplot2::ggplot() +
    geom_jitter(data = df, mapping = aes(x = count, y = total), alpha = .5) +
    stat_ellipse(data = df, mapping = aes(x = count, y = total), type = "norm", alpha = .75,
                 level = data[["chisq"]]) + # Chisq when type = "norm"
    labs(x = "Number of segments",
         y = "Total IBD segment length (cM)") +
    theme_bw()

  if (!is.null(obs.features)) {
    count.obs = obs.features$obs[["count"]]
    total.obs = obs.features$obs[["total"]]

    df.obs = data.frame(count = count.obs, total = total.obs)

    p <- p +
      geom_point(data = df.obs, mapping = aes(x = count, y = total),
                  shape = 4, stroke = 2, size = 4, colour = "red2")
  }

  return (p)

}
