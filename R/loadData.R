#'
#'@export
loadData <- function(data = NULL) {

  if(is.null(data)) {
    dta = prepareData(ibdrel_unilineal)
  } else {
    dta = prepareData(data)
  }

  dta
}

#'
#'@export
prepareData <- function(data) {
  segmentData = data$segments
  peds = data$peds
  metadata = pedsMetadata(peds)

  data.clean = aggregateSegments(segmentData, metadata, "sex.eqclass")

  data.clean
}

#'
#'@export
removePeds <- function(data, peds, mode = "exclude") {
  data.peds <- names(data$peds)

  if (!all(data.peds == names(data$segments))) {
    stop("List of peds and list of segments are not aligned.")
  }

  if(any(!(peds) %in% data.peds)) {
    stop("Relationship is not in data.")
  }

  if (!(mode) %in% c("exclude", "include")) {
    stop("Mode must be either 'exclude' or 'include'.")
  }

  peds.idx <- which(data.peds %in% peds)

  if (mode == "exclude") {
    data$peds <- data$peds[-peds.idx]
    data$segments <- data$segments[-peds.idx]
  } else if (mode == "include") {
    data$peds <- data$peds[peds.idx]
    data$segments <- data$segments[peds.idx]
  }


  data

}
