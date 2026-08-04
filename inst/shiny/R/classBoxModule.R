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

checkSelectedClass = function(selected) {
  if(is.null(selected)) {
    validated(need(FALSE, "Please select a class."))
  }
}

prevPed = function(id) {

  tooltip(
    actionBttn(
      inputId = NS(id, "prevPed"),
      label = NULL,
      size = "sm",
      style = "jelly",
      icon = icon("backward-step")
    ),
    title = "Previous pedigree"
  )
}

nextPed = function(id) {

  tooltip(
    actionBttn(
      inputId = NS(id, "nextPed"),
      label = NULL,
      size = "sm",
      style = "jelly",
      icon = icon("forward-step")
    ),
    title = "Next pedigree"
  )
}

classBoxUI = function(id) {

  bs4Card(
    title = "Plot",
    width = NULL,
    collapse = TRUE,
    collapsed = FALSE,
    status = "olive",

    fluidRow(
      column(
        width = 6,
        uiOutput(NS(id, "pedControls")),
        plotOutput(NS(id, "pedPlot"))
      ),
      column(
        width = 6,
        plotOutput(NS(id,"jointDistPlot"))
      ),
      uiOutput(NS(id, "classSummary"))
    )

    #box(
    #  title = "Pedigrees",
    #  width = NULL,
    #  collapsible = TRUE,
    #  collapsed = TRUE,
    #  div(
    #    class = "ped-plot",
    #    prevPed(id),
    #    nextPed(id)
    #   # plotOutput(NS(id, "pedPlot"))
    #  )
    #),

    #box(
    #  title = "IBD segment distribution",
    #  width = NULL,
    #  collapsible = TRUE,
    #  collapsed = TRUE,
    #  div(
    #    class = "segment-plot"
    #  #  plotOutput(NS(id, "jointDistPlot"))
    #  ),
    #  footer = tagList(
    #    div(
    #      uiOutput(NS(id, "outlierInfo"))
    #    )
    #  )
    # )
  )
}

classBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {


    #observe({
    #  updateSelectInput(session, "selectPedigree",
    #                    choices = subset(data[["metadata"]], class == input$selectClass)[["rel"]])
    #})

    observe({
      req(data[["selectedClass"]])
      #peds = subset(data[["metadata"]], class == data[["selectedClass"]])[["rel"]]
      #updatePickerInput(session, "selectPedigree", choices = peds)

    })

    observeEvent(data[["selectedTab"]], {
      data[["selectedClass"]] = NULL
    })

    # Pedigree plotting

    peds <- reactive({ # Peds in selected class
      metadata <- data[["metadata"]]
      peds_all <- data[["peds"]]
      selected <- subset(metadata, class == data[["selectedClass"]])$rel
      peds_all[selected]
    })

    output$classSummary <- renderUI({
      req(data[["selectedClass"]])

      npeds <- length(peds())

      if (npeds == 1) {
        tags$p("The selected class contains", npeds, "pedigree.")
      } else {
        tags$p("The selected class contains", npeds, "indistinguishable pedigrees.")
      }



    })


    output$outlierInfo <- renderUI({
      req(data[["selectedClass"]])

      checkObs(data[["obs"]])

      makeOutlierOutput(data, data[["selectedTab"]], data[["selectedClass"]])
    })


    output$jointDistPlot = renderPlot({
      req(data[["selectedClass"]])

      resolution = data[["selectedTab"]]

      jointDistPlot(data, resolution, data[["selectedClass"]])
    })





    current_index <- reactiveVal(1)
    observeEvent(c(input$resTabs, data[["selectedClass"]]), {
      current_index(1)
    })


    observeEvent(input$nextPed, {
      idx <- current_index()
      if (idx < length(peds())) current_index(idx + 1)
    })

    output$pedControls <- renderUI({
      req(data[["selectedClass"]])

      tagList(
        prevPed(id),
        helpText(paste0("Pedigree ", current_index(), "/", length(peds()))),
        nextPed(id)
      )
    })


    observeEvent(input$prevPed, {
      idx <- current_index()
      if (idx > 1) current_index(idx - 1)
    })

    output$pedPlot <- renderPlot({

      if (is.null(data[["selectedClass"]])) {
        validate(need(FALSE, "Please select a class from the relatedness table."))
      }


      if (length(peds()) < current_index()) {
        validate(need(FALSE, "Please select a class from the relatedness table."))
      }

      pedPlot(peds()[[current_index()]])

    }, execOnResize = TRUE, res = 72)


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

  mdists = data[["mdists"]][[resolution]]
  thresholds = data[["mdistthreshold"]][[resolution]]

  mdist.idx = which(names(mdists) == selected.class)
  threshold.idx = which(names(thresholds) == selected.class)
  mdist = mdists[mdist.idx]
  threshold = thresholds[threshold.idx]

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
  df <- df[,c("rel", "kappa0", "kappa1", "kappa2", "kinship")]
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

  segments = data[["ct"]][[resolution]][[selected.class]]

  count = segments[["count"]]
  total = segments[["total"]]

  obs = data[["obs"]]

  if (!identical(obs, numeric(0))) {
    obs.features = ibdrel::obsToFeatures(obs, data[["cutoff"]], c("count", "total"), data[["isFeatures"]])
  } else {
    obs.features = NULL
  }
  df = data.frame(count, total)

  p <- ggplot2::ggplot(data = df, mapping = aes(x = count, y = total)) +
    geom_jitter(alpha = .5) +
    stat_ellipse(type = "norm", alpha = .75,
                 level = data[["outlierThreshold"]]) + # Chisq when type = "norm"
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

  #p <- ggExtra::ggMarginal(p, type = "density")

  return (p)

}

# Inspired by `pedbuildr` plot()
pedPlot <- function(ped) {


  ped.leaves = ibdrel::identifyLeaves(ped)
  ped = ped |> pedtools::setSex(leaves, sex= 0)
  verb = verbalisr::verbalise(ped, ids = ped.leaves)[[1]]$rel

  title <- paste0(toupper(substr(verb,1,1)),
                  tolower(substr(verb,2,nchar(verb))))


  plot(ped, hatched = ped.leaves, fill = list(red = ped.leaves),
       labs = NULL, autoScale = TRUE, title = title,
       cex = 1.2, cex.main = 1.5)

}
