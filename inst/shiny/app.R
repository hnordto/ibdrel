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
library(pedtools)

# Load data


segmentDataRel = readRDS(system.file("data", "segments_unilineal_rel.rds", package = "ibdrel"))
segmentDataDon = readRDS(system.file("data", "segments_unilineal_donnelly.rds", package = "ibdrel"))
segmentDataKappa = readRDS(system.file("data", "segments_unilineal_kappa.rds", package = "ibdrel"))
segmentDataKinship = readRDS(system.file("data", "segments_unilineal_kinship.rds", package = "ibdrel"))
segmentDataDeg = readRDS(system.file("data", "segments_unilineal_degree.rds", package = "ibdrel"))

pedsDataRel = readRDS(system.file("data", "peds_unilineal.rds", package = "ibdrel"))
metadata = pedsMetadata(pedsDataRel)

ui <- bs4Dash::bs4DashPage(
  title = "ibdRel",

  header = bs4DashNavbar(
    status = "olive",
    title =
      div(
        class = "apptitle",
        "ibdRel"
      ),
    navbarMenu(
      id = "navmenu",
      navbarTab(tabName = "analysis", text = "Analysis"),
      navbarTab(tabName = "settings", text = "Settings")
    )
  ),



  sidebar = dashboardSidebar(disable = TRUE, minified = FALSE),
  controlbar = NULL,
  footer = NULL,

  body = bs4DashBody(
    useBusyIndicators(spinners = TRUE),
    tabItems(
      tabItem(
        "analysis",
        fluidRow(
          column(
            width = 3,
            h3("Input"),
            textAreaInput("segText", "Segment lengths", rows = 10),
            numericInput("cutoff", "Cutoff", value = 7, min = 0, step = 1),
            radioButtons("classLevel", "Classification level",
                         choices = c("Relationship" = "rel",
                                     "Equivalence class" = "class",
                                     "Kappa" = "kappa",
                                     "Kinship" = "kinship",
                                     "Degree" = "degree"),
                         selected = "rel",
                         inline = TRUE),
            actionButton("classify", "Classify")
          ),
          column(
            width = 9,
            bs4TabCard(
              title = "Analysis",
              collapsible = FALSE,
              width = NULL,
              tabPanel(
                title = "Classification",
                fluidRow(
                  column(
                    width = 6,
                    gt::gt_output("results_table")
                  ),
                  column(
                    width = 6,
                    plotOutput("pedPlot"),
                    uiOutput("choose_ped")
                  )
                )
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
            radioButtons("posteriorProb", "Posterior probability",
                         choices = c("Absolute" = "absolute",
                                     "Relative" = "relative"),
                         selected = "absolute",
                         inline = TRUE),
            radioButtons("normalizedProb", "Probability normalization",
                         choices = c("Normalized" = "normalized",
                                     "Unnormalized" = "unnormalized"),
                         selected = "normalized"),
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
    ),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    )
  )
)
# gt::gt_output("posterior_table")

server <- function(input, output, session) {

  # DATA AND INPUTS -----------------

  trainingData = reactive({
    switch(input$classLevel,
           rel = segmentDataRel,
           class = segmentDataDon,
           kappa = segmentDataKappa,
           kinship = segmentDataKinship,
           degree = segmentDataDeg)
  })

  aggLevel = reactive({input$classLevel})

  metadata = reactive({
    pedsMetadata(pedsDataRel)
  })

  peds = reactive({
    switch(input$classLevel,
           rel = pedsDataRel,
           class = NULL,
           kappa = NULL,
           kinship = NULL,
           degree = NULL)
  })

  features = reactive({
    lapply(trainingData(), prepareFeatures, cutoff = input$cutoff)
  })

  pdfuns = reactive({
    req(trainingData())
    lapply(trainingData(), preparePdfs, cutoff = input$cutoff)
  })
  # Force this to run immediately
  observe({
    pdfuns()
  })

  obs = reactive({
    as.numeric(input$segText |> strsplit("\n") |> unlist())
  })

  output$choose_ped_var <- renderUI({
    selectInput("choose_ped_var", "Class", choices = names(posteriors()))
  })

  output$choose_ped <- renderUI({
    selectInput("choose_ped", "Choose pedigree", choices = names(posteriors()))
  })

  # CALCULATIONS -----------------

  posteriors = reactiveVal(NULL)
  lofs = reactiveVal(NULL)
  outliers = reactiveVal(NULL)
  observeEvent(input$classify, {
    post = classify(obs(), pdfuns())
    posteriors(post)

    lof = LOF(obs(), features())
    lofs(lof$lof)

    outlier = ifelse(lof$lof > lof$threshold, TRUE, FALSE)
    outliers(outlier)

  })

  # OUTPUTS -----------------

  # Plots
  output$varScatterplot <- renderPlot({

    req(input$name_clicked)

    outlierPlot(features(),
                obsToFeatures(obs()),
                input$name_clicked)
  })

  output$pedPlot <- renderPlot({

      validate(
        need(!is.null(peds()), "Pedigree plotting not supported for chosen settings.")
      )

      req(peds, input$choose_ped)

      ped = peds()[input$choose_ped]
      ped.leaves = identifyLeaves(ped)

      tryCatch(
        plot(ped, cex = input$cex),
        error = function(e) {
          plot.new()
          msg = if(grepl("reduce cex", conditionMessage(e))) "(Too big for plot region. Reduce cex in 'Settings'.)" else conditionMessage(e)
          mtext(msg, line = 0, col = 2)
        }
      )

      plot(ped, cex = input$cex,
           labs = NULL,
           hatched = ped.leaves,
           fill = list(red = ped.leaves))
    })



  # Tables

  results_tbl = eventReactive(input$classify, {



    tbl <- resultTable(metadata(),
                       posteriors(),
                       outliers(),
                       lofs(),
                       length(pdfuns()[[1]])-1,
                       aggLevel())
    tbl
  })

  output$eqclass_list = renderPrint({
    req(input$name_clicked)


  })

  output$rels_list = renderPrint({
    req(input$name_clicked)

    relsAtLevel(metadata(),
                input$classLevel,
                input$name_clicked)

  })

  output$results_table = render_gt({
    req(results_tbl())
  })
  output$results_table_full = render_gt({
    req(results_tbl_full())
  })

  # Modals
  observeEvent(input$name_clicked, {
    showModal(modalDialog(
      title = input$name_clicked,
      h5("Relationships in class"),
      uiOutput("rels_list"),
      h5("Training data distribution (2D)"),
      uiOutput("choose_var1"),
      uiOutput("choose_var2"),
      plotOutput("varScatterplot")
    ))
  })

  output$test = renderPrint(posteriors())


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
