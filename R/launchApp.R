#' Launch ibdRel
#'
#' This launchaes the ibdRel shiny applications
#'
#' @export

launchApp = function() {
  suppressPackageStartupMessages(
    shiny::runApp(system.file("shiny", package = "ibdrel"), launch.browser = TRUE)
  )
}
