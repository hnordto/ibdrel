#'
#'@export
outlierAnalysis <- function(object, ...) {
  UseMethod("outlierAnalysis")
}

#'
#'@export
outlierAnalysis.ibdrelModel <- function(object, newdata, threshold = 0.95, collapseDistant = TRUE) {

  segments = object$segments
  metadata = object$metadata
  cutoff = object$cutoff

  features.count.total <- extractCountTotal(model)

  if(!is.null(names(get_inner_element(newdata)))) {
    if(!any(names(get_inner_element(newdata)) %in% featureSel)) {
      stop("Input summaries do not match model features.")
    }

    outlier <- lapply(features.count.total, function(x) {
      lapply(newdata, function(obs) {
        distance(obs = obs, features = x, cutoff = cutoff, featureSel = featureSel, threshold = threshold, isFeatures = TRUE)
      })
    })

  } else {
    allobs = get_inner_element(newdata)

    outlier <- lapply(features.count.total, function(x) {
      lapply(allobs, function(obs) {
        distance(obs = obs, features = x, cutoff = cutoff, featureSel = featureSel, threshold = threshold, isFeatures = FALSE)
      })
    })
  }

  return(outlier)

}

#'
#'@export
extractCountTotal <- function(object, ...) {
  UseMethod("extractCountTotal")
}

#'
#'@export
extractCountTotal.ibdrelModel <- function(object, ...) {
  segments = object$segments
  metadata = object$metadata
  cutoff = object$cutoff

  featureSel = c("count", "total")

  segments.agg <- lapply(RESOLUTIONS, function(x) {
    aggregateSegments(segments, metadata, x, collapseDistant)
  })
  names(segments.agg) <- RESOLUTIONS

  features.count.total <- lapply(segments.agg, function(x) {
    lapply(x, prepareFeatures, featureSel = featureSel, cutoff = cutoff)
  })

  return (features.count.total)
}

reportOutlier.ibdrelModel <- function(object, newdata, threshold = 0.95, collapseDistant = TRUE) {

  outlier <- outlierAnalysis.ibdrelModel(object, newdata, threshold, collapseDistant)

  report <- lapply(outlier, function(res) {
    lapply(res, function(obs) {
      lapply(obs, function(rel) {

        if(is.na(rel$dist) || is.na(rel$threshold)) {
          NA
        } else {
          if (rel$dist > rel$threshold) TRUE else FALSE
        }
      })
    })
  })

  return (report)

}
