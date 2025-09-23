helpWidget = function(id) {
  tooltip(
    actionBttn(
      inputId = NS(id, "help"),
      label = NULL,
      icon = icon("question")
    ),
    title = "See help on this panel."
  )
}

checkLikelihood = function(likelihood) {
  if (all(likelihood == -Inf)) {
    validate(need(FALSE, "Waiting for input.\nNeed at least one segment with length >= cutoff."))
  }
}

likelihoodBoxUI = function(id) {


  tabBox(
    id = NS(id, "resTabs"),
    width = NULL,
    selected = "eqclass.detailed",
    collapsible = FALSE,

    title = div(
      div(
        "Result table",
      ),
      div(
        helpWidget(id)
      )
    ),

    tabPanel(
      title = "Sex-specific pedigree class",
      value = "eqclass.detailed",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableDetailedEq"))
      )
    ),

    tabPanel(
      title = "Pedigree class",
      value = "eqclass",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableEq"))
      )
    ),

    tabPanel(
      title = "Kappa coefficients",
      value = "kappa",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableKappa"))
      )
    ),

    tabPanel(
      title = "Kinship coefficient",
      value = "kinship",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableKinship"))
      )
    ),

    tabPanel(
      title = "Degree",
      value = "degree",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableDegree"))
      )
    )
  )
}

likelihoodBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    # Set active tab
    observe({
      req(input$resTabs)
      data[["selectedTab"]] <- input$resTabs
    })

    output$resTableDetailedEq = gt::render_gt({
      res <- "eqclass.detailed"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res)
    })

    output$resTableEq = gt::render_gt({
      res <- "eqclass"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res)
    })

    output$resTableKappa = gt::render_gt({
      res <- "kappa"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res)
    })

    output$resTableKinship = gt::render_gt({
      res <- "kinship"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res)
    })

    output$resTableDegree = gt::render_gt({
      res <- "degree"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res)
    })


    # Help
    observeEvent(input$help, {
      shinyalert(
        className = "helpbox",
        html = TRUE,
        text = read_file("help/likelihoodBox.html"),
        showConfirmButton = FALSE,
        closeOnClickOutside = TRUE
      )
    })

  })
}

# Format table

likelihoodTable <- function(data, resolution) {

  likelihoods = data[["likelihoods"]][[resolution]] # XX FIX XX
  mdists = data[["mdists"]][[resolution]]
  filter = data[["filter"]]
  normalize = data[["normalize"]]
  threshold = data[["mdistthreshold"]]


  if (filter) {
    inlier.classes <- names(which(mdists < threshold))
    likelihoods <- likelihoods[inlier.classes]
  }

  if (normalize) {
    likelihoods = ibdrel::normalizeLikelihoods(likelihoods)
  }

  likelihoods <- sort(likelihoods, decreasing = TRUE)
  classlabels <- sapply(names(likelihoods), ibdrel::classTranslator, resolution)


  df <- data.frame(Relationship = classlabels,
                   Likelihood = round(likelihoods, 2))
  df <- gt(df) |>
    gt_theme_538(quiet = TRUE)
  return (df)
}
