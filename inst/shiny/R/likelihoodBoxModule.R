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

checkInput = function(data) {
  if (is.null(data)) {
    validate(need(FALSE, "Waiting for input."))
  }
  data
}

likelihoodBoxUI = function(id) {


  tabBox(
    id = NS(id, "resTabs"),
    width = NULL,
    selected = "eqclass.detailed",

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

      likelihoods <- checkInput(data[["likelihoods"]])

      input = if (data[["normalize"]]) {
        likelihoods[["norm"]][["eqclass.detailed"]]
      } else {
        likelihoods[["raw"]][["eqclass.detailed"]]
      }

      names.input = sapply(names(input), ibdrel::classTranslator, "eqclass.detailed")

      likelihoodTable(names.input, input)
    })

    output$resTableEq = gt::render_gt({

      likelihoods <- checkInput(data[["likelihoods"]])

      input = if (data[["normalize"]]) {
       likelihoods[["norm"]][["eqclass"]]
      } else {
        likelihoods[["raw"]][["eqclass"]]
      }

      names.input = sapply(names(input), ibdrel::classTranslator, "eqclass")

      likelihoodTable(names.input, input)
    })

    output$resTableKappa = gt::render_gt({

      likelihoods <- checkInput(data[["likelihoods"]])

      input = if (data[["normalize"]]) {
        likelihoods[["norm"]][["kappa"]]
      } else {
        likelihoods[["raw"]][["kappa"]]
      }

      names.input = sapply(names(input), ibdrel::classTranslator, "kappa")

      likelihoodTable(names.input, input)
    })

    output$resTableKinship = gt::render_gt({

      likelihoods <- checkInput(data[["likelihoods"]])

      input = if (data[["normalize"]]) {
        likelihoods[["norm"]][["kinship"]]
      } else {
        likelihoods[["raw"]][["kinship"]]
      }

      names.input = sapply(names(input), ibdrel::classTranslator, "kinship")

      likelihoodTable(names.input, input)
    })

    output$resTableDegree = gt::render_gt({

      likelihoods <- checkInput(data[["likelihoods"]])

      input = if (data[["normalize"]]) {
        likelihoods[["norm"]][["degree"]]
      } else {
        likelihoods[["raw"]][["degree"]]
      }

      names.input = sapply(names(input), ibdrel::classTranslator, "degree")

      likelihoodTable(names.input, input)
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

likelihoodTable <- function(rels, likelihoods) {
  df <- data.frame(Relationship = rels,
                   Likelihood = round(likelihoods, 2))
  df <- gt(df) |>
    gt_theme_538()
  return (df)
}
