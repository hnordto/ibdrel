library(shiny)
library(bs4Dash)
library(ggplot2)
library(stringr)
library(gt)
library(DT)
library(scales)
library(dplyr)
library(gtExtras)
library(pedtools)

# Load data


segmentDataRel = readRDS(system.file("data", "segments_unilineal.rds", package = "ibdrel"))
segmentDataDeg = readRDS(system.file("data", "segments_unilineal_deg.rds", package = "ibdrel"))

pedsDataRel = readRDS(system.file("data", "peds_unilineal.rds", package = "ibdrel"))
pedsMetadataRel = readRDS(system.file("data", "peds_metadata_unilineal.rds", package = "ibdrel"))


ui <- bs4Dash::bs4DashPage(
  title = "ibdClassifier",

  header = bs4DashNavbar(
    status = "olive",
    title = "ibdClassifier",
    navbarMenu(
      id = "navmenu",
      navbarTab(tabName = "analysis", text = "Analysis"),
      navbarTab(tabName = "settings", text = "Settings")
    )
  ),



  sidebar = dashboardSidebar(disable = TRUE, minified = FALSE),

  body = bs4DashBody(
    tabItems(
      tabItem(
        "analysis",
        fluidRow(
          column(
            width = 4,
            h3("Input"),
            textAreaInput("segText", "Segment lengths", rows = 10),
            numericInput("cutoff", "Cutoff", value = 0, min = 0, step = 1),
            numericInput("num_results", "Number of results to display",
                         value = 5, min = 1),
            actionButton("classify", "Classify")
          ),
          column(
            width = 8,
            bs4TabCard(
              title = "Analysis",
              collapsible = FALSE,
              width = NULL,
              tabPanel(
                title = "Overview",
                fluidRow(
                  column(
                    width = 6,
                    gt::gt_output("results_table"),
                  ),
                  column(
                    width = 6,
                    plotOutput("pedPlot"),
                    uiOutput("choose_ped")
                  )
                )
              ),
              tabPanel(
                title = "Variable plots",
                uiOutput("choose_ped_var"),
                uiOutput("choose_var1"),
                uiOutput("choose_var2"),
                plotOutput("varScatterplot")
              ),
              tabPanel(
                title = "Full posterior table",
                gt::gt_output("results_table_full")
              )
            )
          )
        )
      ),
      tabItem(
        "settings",
        fluidRow(
          column(
            width = 4,
            h3("Settings"),
            radioButtons("classLevel", "Classification level",
                         choices = c("Automatic" = "auto",
                                     "Degree" = "degree",
                                     "Relationship" = "relationship",
                                     "Sex-specific relationship" = "sexspecific"),
                         selected = "auto",
                         inline = TRUE),
            radioButtons("posteriorProb", "Posterior probability",
                         choices = c("Absolute" = "absolute",
                                     "Relative" = "relative"),
                         selected = "absolute",
                         inline = TRUE),
            numericInput("outlier_threshold", "Chi-square Mahalanobis outlier threshold",
                         value = 0.05),
            numericInput("cex", "CEX", value = 1)
          ),
          column(
            width = 8,
            bs4TabCard(
              title = "Model",
              collapsible = FALSE,
              width = NULL,
              tabPanel(
                title = "Model",
                DTOutput("class_tbl"),
                actionButton("update_classes", "Update class selection")
              )
            )
          )
        )
      )
    )
  )
)
# gt::gt_output("posterior_table")

server <- function(input, output, session) {

  # DATA AND INPUTS -----------------

  trainingData = reactive({
    switch(input$classLevel,
           auto = segmentDataRel,
           degree = segmentDataDeg,
           relationship = segmentDataRel) # not supported yet
  })

  peds = reactive({
    switch(input$classLevel,
           auto = pedsDataRel,
           degree = NULL,
           relationship = pedsDataRel)
  })

  pedsMetadata = reactive({
    switch(input$classLevel,
           auto = pedsMetadataRel,
           degree = NULL,
           relationship = pedsMetadataRel)
  })

  features = reactive({
    lapply(trainingData(), prepareFeatures)
  })

  pdfuns = reactive({
    lapply(trainingData(), preparePdfs)
  })

  obs = reactive({
    as.numeric(input$segText |> strsplit("\n") |> unlist())
  })

  output$choose_ped_var <- renderUI({
    selectInput("choose_ped_var", "Class", choices = names(posteriors()))
  })

  output$choose_var1 <- renderUI({
    selectInput("choose_var1", "Variable 1", choices = names(pdfuns()[[1]]))
  })

  output$choose_var2 <- renderUI({
    selectInput("choose_var2", "Variable 2", choices = names(pdfuns()[[1]]))
  })

  output$choose_ped <- renderUI({
    selectInput("choose_ped", "Choose pedigree", choices = names(posteriors()))
  })

  # CALCULATIONS -----------------

  posteriors = reactiveVal(NULL)
  mdists = reactiveVal(NULL)
  outliers = reactiveVal(NULL)
  observeEvent(input$classify, {
    post = classify(obs(), pdfuns())
    posteriors(post)

    mdist = distance(obs(), features())
    mdists(mdist)

    outlier = ifelse(mdist > qchisq(p = 1-input$outlier_threshold,
                                    df = length(pdfuns()[[1]])-1), T, F)
    outliers(outlier)

  })

  # OUTPUTS -----------------

  # Plots
  output$varScatterplot <- renderPlot({

    req(input$choose_ped_var, input$choose_var1, input$choose_var2)

    varScatterplot(features(),
                   obsToFeatures(obs()),
                   input$choose_ped_var,
                   input$choose_var1,
                   input$choose_var2)
  })

  output$pedPlot <-renderPlot({

      validate(
        need(!is.null(peds()), "Pedigree plotting not supported for chosen settings.")
      )

      req(peds, input$choose_ped)

      ped = peds()[input$choose_ped]

      tryCatch(
        plot(ped, cex = input$cex),
        error = function(e) {
          plot.new()
          msg = if(grepl("reduce cex", conditionMessage(e))) "(Too big for plot region. Reduce cex in 'Settings'.)" else conditionMessage(e)
          mtext(msg, line = 0, col = 2)
        }
      )

      plot(ped, cex = input$cex)
    })


  # Tables

  results_tbl = eventReactive(input$classify, {
    tbl <- resultTable(pedsMetadata(),
                       posteriors(),
                       outliers(),
                       mdists(),
                       length(pdfuns()[[1]])-1,
                       input$num_results)
    tbl |> cols_hide(columns = c(Distance, Distance_p, Group, class))
  })

  results_tbl_full = eventReactive(input$classify, {
    tbl <- resultTable(pedsMetadata(),
                       posteriors(),
                       outliers(),
                       mdists(),
                       length(pdfuns()[[1]])-1,
                       NULL)
    tbl |> cols_hide(columns = c(Group, class))
  })

  output$results_table = render_gt({
    req(results_tbl())
  })
  output$results_table_full = render_gt({
    req(results_tbl_full())
  })

  # SETINGS -----------------

  class_df <- reactive({
    data.frame(
      class = names(trainingData()),
      check = shinyInput(checkboxInput, length(trainingData()), "checkb")
    )
  })
  output$class_tbl <- renderDT({
    datatable(
      class_df(),
      rownames = FALSE,
      escape = FALSE,
      callback = JS(c(
        "$('[id^=checkb]').on('click', function(){",
        "  var id = this.getAttribute('id');",
        "  var i = parseInt(/checkb(\\d+)/.exec(id)[1]);",
        "  var value = $(this).prop('checked');",
        "  var info = [{row: i, col: 3, value: value}];",
        "  Shiny.setInputValue('dtable_cell_edit:DT.cellInfo', info);",
        "})"))
    )
  })

}

shinyApp(ui, server)
