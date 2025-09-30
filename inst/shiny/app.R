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

    tags$script(HTML("
    Shiny.addCustomMessageHandler('initClickableCells', function(message) {
  var inputId = message.ns; // namespaced inputId
  $('.clickable-cell').off('click').on('click', function() {
    var val = $(this).attr('id').replace(/^rel_/, '');
    Shiny.setInputValue(inputId, val, {priority:'event'});

    // Highlight clicked cell
    $('.clickable-cell').css('background-color', '');
    $(this).css('background-color', '#ffeaa7');
  });
});;
  ")),


    fluidRow(

      # Input
      column(
        width = 3,
        inputBoxUI("input")
      ),

      # Likelihood table
      column(
        width = 4,

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
