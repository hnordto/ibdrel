#'
#'@export
loadData <- function(data = "ibdrel-unilineal",
                     levels = "all") {
  if (length(levels) == 1 && levels == "all") {
    levels = RESOLUTIONS
  } else if (any(!(levels %in% RESOLUTIONS))) {
    stop("Invalid level. Must be one or multiple of ",
         paste(RESOLUTIONS, collapse = " , "))
  }

  if(!(data) %in% BULTIN.MODELS) {
    stop("Unknown model: ", model)
  }

  if (data == "ibdrel-unilineal") {
    segmentData = ibdrel_unilineal$segments
    peds = ibdrel_unilineal$peds
    metadata = pedsMetadata(peds)
  }

  data = sapply(RESOLUTIONS, function(x) {
    aggregateSegments(segmentData, metadata, x)
  }, simplify = FALSE)

  return (data)
}
