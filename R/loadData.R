#'
#'@export
loadData <- function(data = NULL,
                     levels = "all",
                     conditionDistant = FALSE) {
  if (length(levels) == 1 && levels == "all") {
    levels = RESOLUTIONS
  } else if (any(!(levels %in% RESOLUTIONS))) {
    stop("Invalid level. Must be one or multiple of ",
         paste(RESOLUTIONS, collapse = " , "))
  }

  if(is.null(data)) {
    dta = prepareData(ibdrel_unilineal, levels, conditionDistant)
  } else {
    dta = prepareData(data, levels, conditionDistant)
  }

  dta
}

#'
#'@export
prepareData <- function(data, levels = "all", conditionDistant = F) {
  segmentData = data$segments
  peds = data$peds
  metadata = pedsMetadata(peds)

  data.clean = sapply(RESOLUTIONS, function(x) {
    aggregateSegments(segmentData, metadata, x, conditionDistant)
  }, simplify = FALSE)

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
