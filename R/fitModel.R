
fitModel <- function(data,
                     features = c("count", "total"),
                     cutoff = 7) {
  if (is.null(data)) {
    data = loadData()
  }

  if (any(!(features %in% VALID.FEATURES))) {
    stop("Invalid feature(s). Must be one or multiple of ",
         paste(VALID.FEATURES, collapse = ", "))
  }

  pdfs = lapply(data, function(x) {
    lapply(x, preparePdfs,
           featureSel = features,
           cutoff = cutoff)
  })

  model = list(
    pdfs = pdfs,
    cutoff = cutoff,
    features = features
  )

  return (model)

}
