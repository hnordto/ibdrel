likelihoodBoxUI = function(id) {

  tabBox(
    id = NS(id, "resTabs"),
    width = NULL,
    selected = "Sex-specific pedigree class",

    title = "Result table",

    tabPanel(
      title = "Sex-specific pedigree class",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableDetailedEq"))
      )
    ),

    tabPanel(
      title = "Pedigree class",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableEq"))
      )
    ),

    tabPanel(
      title = "Kappa coefficients",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableKappa"))
      )
    ),

    tabPanel(
      title = "Kinship coefficient",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableKinship"))
      )
    ),

    tabPanel(
      title = "Degree",
      div(
        class = "table-box",
        gt::gt_output(NS(id, "resTableDegree"))
      )
    ),

    materialSwitch(inputId = NS(id, "normalizeButton"), value = TRUE,
                   label = "Normalize likelihoods?"),
    materialSwitch(inputId = NS(id, "filterButton"), value = FALSE,
                   label = "Filter outliers?")
  )
}

likelihoodBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    output$resTableDetailedEq = gt::render_gt({
      input = if (input$normalizeButton) {
        data[["likelihoods"]][["norm"]][["eqclass.detailed"]]
      } else {
        data[["likelihoods"]][["raw"]][["eqclass.detailed"]]
      }

      likelihoodTable(names(input), input)
    })

    output$resTableEq = gt::render_gt({
      input = if (input$normalizeButton) {
        data[["likelihoods"]][["norm"]][["eqclass"]]
      } else {
        data[["likelihoods"]][["raw"]][["eqclass"]]
      }

      likelihoodTable(names(input), input)
    })

    output$resTableKappa = gt::render_gt({
      input = if (input$normalizeButton) {
        data[["likelihoods"]][["norm"]][["kappa"]]
      } else {
        data[["likelihoods"]][["raw"]][["kappa"]]
      }

      likelihoodTable(names(input), input)
    })

    output$resTableKinship = gt::render_gt({
      input = if (input$normalizeButton) {
        data[["likelihoods"]][["norm"]][["kinship"]]
      } else {
        data[["likelihoods"]][["raw"]][["kinship"]]
      }

      likelihoodTable(names(input), input)
    })

    output$resTableDegree = gt::render_gt({
      input = if (input$normalizeButton) {
        data[["likelihoods"]][["norm"]][["degree"]]
      } else {
        data[["likelihoods"]][["raw"]][["degree"]]
      }

      likelihoodTable(names(input), input)
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
