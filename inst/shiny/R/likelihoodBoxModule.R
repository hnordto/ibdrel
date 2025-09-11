likelihoodBoxUI = function(id) {
  box(
    id = NS(id, "box"),
    width = NULL,
    title = "Likelihood table",
    div(
      class = "table-box",
      gt_output(NS(id, "likelihoodTable"))
    ),
    materialSwitch(inputId = NS(id, "normalizeButton"), value = TRUE,
                label = "Normalize likelihood?")
  )
}

likelihoodBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    output$likelihoodTable = render_gt({
      req(data[["likelihoods"]])

      likelihoods = if (input$normalizeButton) {
        data[["likelihoods"]][["norm"]]
      } else {
        data[["likelihoods"]][["raw"]]
      }

      likelihoods = likelihoods
      rels.ordered = names(likelihoods)
      likelihoodTable(rels.ordered,
                      likelihoods)
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
