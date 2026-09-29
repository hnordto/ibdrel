#'
#' @export
predict.ibdrelModel <- function(object, newdata,
                    normalize = TRUE, sort = TRUE, collapseDistant = TRUE) {

  pdfs = object$pdfs

  allobs = get_inner_element(newdata)

  if (length(obs) > 1 && isTRUE(sort)) {
    stop("Cannot sort multiple observations, as ranks may differ. Set sort = FALSE.")
  }

  pred = lapply(allobs, function(obs) {
    classify(obs, pdfuns = pdfs, cutoff = object$cutoff, sort = sort)
  })
  pred = do.call(rbind, pred)

  pred <- lapply(RESOLUTIONS, function(x) {
    aggregateProbs(pred, pedsMetadata(model$peds), x, collapseDistant)
  })
  names(pred) <- RESOLUTIONS


  if (normalize) {
    pred = lapply(pred, function(mat) {
      t(apply(mat, 1, normalizeLikelihoods))
    })
  }


  return (pred)
}

aggregateProbs <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {
  lookup = lookupClass(NULL, metadata.agg.column, "rel", metadata, collapseDistant)

  if (isTRUE(collapseDistant)) {
    lookup["distant"] = "distant"
  }

  agg.ids = lookup[colnames(pred)]

  t(rowsum(t(pred), group = agg.ids))
}

# Helper function to retrieve each observation when multiple are inputted
get_inner_element <- function(x) {
  if (!is.list(x)) {
    return(list(x))
  } else {
    return(unlist(lapply(x, get_inner_element), recursive = FALSE))
  }
}
