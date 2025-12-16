#' Launch the ibdrel app
#'
#' This launches the Shiny application for performing pairwise relationship estimation.
#'
#' @examples
#' \dontrun{
#' launchApp()
#' }
#'
#' @export

launchApp = function() {
  suppressPackageStartupMessages(
    shiny::runApp(system.file("shiny", package = "ibdrel"), launch.browser = TRUE)
  )
}
