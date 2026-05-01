ui = dashboardPage(

  # Header
  header = dashboardHeader(

    status = "olive",

    title =
      div(
        class = "apptitle",
        "ibdrel"
      )
  ),

  # Sidebar
  sidebar = dashboardSidebar(disable = TRUE, minified = FALSE, width = 0),

  # Body

  body = dashboardBody(

    useShinyjs(),
    useBusyIndicators(),

    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),


    fluidRow(

      # Input
      column(
        width = 2,
        inputBoxUI("input")
      ),

      # Likelihood table
      column(
        width = 5,

        likelihoodBoxUI("likelihoodTable")
      ),

      # Class analysis
      column(
        width = 5,
        classBoxUI("classInfo")
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

}

shinyApp(ui, server)
