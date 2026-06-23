ui = bs4DashPage(

  dark = NULL,
  help = NULL,

  # Header
  header = bs4DashNavbar(

    status = "olive",

    title =
      div(
        class = "apptitle",
        "IBDrel"
      ),

    rightUi =  tagList(tags$li(class = "nav-item dropdown",
      div(class = "aligned-row", style = "margin-right: 22.5px; gap: 15px;",
          actionBttn("settings", icon("gear"), style = "jelly", color = "danger", size = "m")
      )
    ))
  ),

  # Sidebar
  sidebar = bs4DashSidebar(disable = TRUE, minified = FALSE, width = 0),

  # Body

  body = bs4DashBody(

    useShinyjs(),
    useBusyIndicators(),

    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),

    fluidRow(
      column(
        width = 6,
        inputBoxUI("input"),
        classBoxUI("classInfo")
      ),
      column(
        width = 6,
        likelihoodBoxUI("likelihoodTable")
      )
    )
  )
)

server <- function(input, output, session) {

  data = reactiveValues()
  data[["obs"]] = NULL

  inputBoxServer("input", data)
  likelihoodBoxServer("likelihoodTable", data)
  classBoxServer("classInfo", data)

  observeEvent(input$loadfile, {
    x = scan(input$loadfile$datapath)
  })

  settings = reactiveValues(cutoff = 7,
                            inputType = "segments",
                            featureSelection = c("count", "total"),
                            outlierThreshold = 0.95)

  # Settings
  observeEvent(input$settings, {
    showModal(modalDialog(
      div(
        h3("Settings"),
        h5("Match threshold"),
        numericInput("cutoff", NULL,
                     value = settings$cutoff, min = 0, step = 1, width = "auto"),
        hr(),
        h5("Input type"),
        awesomeRadio("inputType", NULL,
                     choices = input.types, selected = settings$inputType,
                     inline = TRUE, width = "auto"),
        hr(),
        h5("Outlier detection"),
        numericInput("outlierThreshold", NULL,
                     value = settings$outlierThreshold, min = 0, max = 1, step = 0.01),
        hr(),
        h5("Feature selection"),
        pickerInput("featureSelection", NULL,
                    choices = c("count", "total"),
                    selected = c("count", "total"),
                    options = pickerOptions(
                      actionsBox = TRUE,
                      size = 10,
                      selectedTextFormat = "count > 3"
                    ), multiple = TRUE)

      ),
      easyClose = TRUE,
      footer = modalButton("Save and close")
    ))
  })

  observeEvent(input$cutoff, {settings$cutoff = input$cutoff})
  observeEvent(input$inputType, {settings$inputType = input$inputType})
  observeEvent(input$featureSelection, {settings$featureSelection = input$featureSelection})
  observeEvent(input$outlierThreshold, {settings$outlierThreshold = input$outlierThreshold})

  observe({
    data[["cutoff"]] = settings$cutoff
    data[["inputType"]] = settings$inputType
    data[["featureSelection"]] = settings$featureSelection
    data[["isFeatures"]] = ifelse(settings$inputType == "segments", F, T)
    data[["outlierThreshold"]] = settings$outlierThreshold
  })
}

shinyApp(ui, server)
