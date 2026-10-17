#'
#' @export
predict.ibdrelModel <- function(object, newdata,
                    metadata = NULL, normalize = TRUE, sort = TRUE, collapseDistant = TRUE) {

  # Pre-loaded metadata speeds up prediction

  pdfs = object$pdfs

  allobs = get_inner_element(newdata)

  if (length(allobs) > 1 && isTRUE(sort)) {
    stop("Cannot sort multiple observations, as ranks may differ. Set sort = FALSE.")
  }

  pred = lapply(allobs, function(obs) {
    classify(obs, pdfuns = pdfs, cutoff = object$cutoff, sort = sort)
  })
  pred = do.call(rbind, pred)

  if (is.null(metadata)) {
    metadata = pedsMetadata(object$peds)
  }


  if (normalize) {
    pred = t(apply(pred, 1, normalizeLikelihoods))


    pred <- lapply(RESOLUTIONS, function(x) {
      aggregateProbs(pred, metadata, x, collapseDistant)
    })
    names(pred) <- RESOLUTIONS
  } else {
    pred <- lapply(RESOLUTIONS, function(x) {
      aggregateLogs(pred, metadata, x, collapseDistant)
    })
    names(pred) <- RESOLUTIONS
  }

  return (pred)
}

aggregateLookup <- function(metadata, metadata.agg.column, collapseDistant) {
  lookup = lookupClass(NULL, metadata.agg.column, "rel", metadata, collapseDistant)

  if (isTRUE(collapseDistant)) {
    lookup["distant"] = "distant"
  }

  agg.ids = lookup[colnames(pred)]

  agg.ids
}

aggregateProbs <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {

  agg.ids = aggregateLookup(metadata, metadata.agg.column, collapseDistant)

  t(rowsum(t(pred), group = agg.ids))
}

# Perhaps numerically unstable
aggregateLogs <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {
  agg.ids = aggregateLookup(metadata, metadata.agg.column, collapseDistant)

  t(log(rowsum(t(exp(pred)), group = agg.ids)))
}

# Helper function to retrieve each observation when multiple are inputted
get_inner_element <- function(x) {
  if (!is.list(x)) {
    return(list(x))
  } else {
    return(unlist(lapply(x, get_inner_element), recursive = FALSE))
  }
}
