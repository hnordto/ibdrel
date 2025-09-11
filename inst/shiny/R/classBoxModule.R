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
      tags$h5("Outlier diagnostics")
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

    observe({
      updatePickerInput(session, "selectClass", choices = data[["classes"]])
      #data[["selectedClass"]] = input$selectClass
    })

    observe({
      updateSelectInput(session, "selectPedigree",
                        choices = subset(data[["metadata"]], class == input$selectClass)[["rel"]])
    })

    output$classInfo <- renderUI({
      req(input$selectClass)

      if (is.null(data[["likelihoods"]])) {
        return (
          div(class = "no-data",
            "Waiting for data."
          )
        )
      }

      makeSummaryOutput(data[["metadata"]],
                        data[["likelihoods"]],
                        input$selectClass)
    })

    output$relationshipTable = render_gt({
      relationshipTable(data[["metadata"]], input$selectClass)
    })

    output$jointDistPlot = renderPlot({
      if (!is.null(data[["obs"]])) {
        jointDistPlot(train.features = data[["features"]][[input$selectClass]],
                      class = input$selectClass,
                      obs.features = ibdrel::obsToFeatures(data[["obs"]]))
      } else {
        jointDistPlot(train.features = data[["features"]][[input$selectClass]],
                      class = input$selectClass,
                      obs.features = NULL)
      }
    })


  })
}

makeSummaryOutput <- function(metadata, likelihoods, selected.class) {

  met = subset(metadata, class == selected.class)

  tagList(
    tags$p(
      tags$b("Selected class:"), tags$i(ibdrel::classTranslator(selected.class)), "\n",
      tags$ul(
        tags$li(tags$b("Log-likelihood:"), likelihoods[["raw"]][[selected.class]]),
        tags$li(tags$b("Normalized likelihood:"), likelihoods[["norm"]][[selected.class]]),
        tags$li(tags$b("Rank:"), which(names(likelihoods[["raw"]]) == selected.class)),
        tags$li(tags$b("Number of pedigrees in class:"), length(met[["code"]]))
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
    gt_theme_538() |>
    cols_label(
      kappa0 = html("&kappa;<sub>0</sub>"),
      kappa1 = html("&kappa;<sub>1</sub>"),
      kappa2 = html("&kappa;<sub>2</sub>"),
      kinship = html("&phi;")
    )
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
