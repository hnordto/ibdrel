#'
#' @export
predict.ibdrelModel <- function(object, newdata,
                    metadata = NULL, normalize = TRUE, sort = TRUE, collapseDistant = TRUE,
                    inputType = "segments") {

  # Pre-loaded metadata speeds up prediction

  pdfs = object$pdfs
  cutoff = object$cutoff
  metadata = if (is.null(metadata)) object$metadata

  # Recognize input type: Vector/unnamed list (segments), named list (summaries)

  if(inputType == "summaries") {

    if(!all(names(get_inner_element(newdata)) %in% object$features)) {
      stop("Input summaries do not match model features.")
    }

    pred = lapply(newdata, function(obs) {
      classify(obs, pdfuns = pdfs, cutoff = cutoff, sort = sort, obsType = "summaries")
    })

  } else {
    allobs = get_inner_element(newdata)

    if (length(allobs) > 1 && isTRUE(sort)) {
      stop("Cannot sort multiple observations, as ranks may differ. Set sort = FALSE.")
    }

    pred = lapply(allobs, function(obs) {
      classify(obs, pdfuns = pdfs, cutoff = cutoff, sort = sort)
    })
  }

  pred = do.call(rbind, pred)


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

aggregateLookup <- function(pred, metadata, metadata.agg.column, collapseDistant) {

  lookup = lookupClass(NULL, metadata.agg.column, "rel", metadata, collapseDistant)

  if (isTRUE(collapseDistant)) {
    lookup["distant"] = "distant"
  }

  agg.ids = lookup[colnames(pred)]

  agg.ids
}

aggregateProbs <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {

  agg.ids = aggregateLookup(pred, metadata, metadata.agg.column, collapseDistant)

  t(rowsum(t(pred), group = agg.ids))
}

# Perhaps numerically unstable
aggregateLogs <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {
  agg.ids = aggregateLookup(pred, metadata, metadata.agg.column, collapseDistant)

  t(log(rowsum(t(exp(pred)), group = agg.ids)))
}

# Helper function to retrieve each observation when multiple are inputted

get_inner_element <- function(x) {
  if (!is.list(x)) {
    return(list(x))
  }

  result <- list()

  for (i in seq_along(x)) {
    element <- x[[i]]
    inner <- get_inner_element(element)

    if (!is.list(element) && !is.null(names(x))) {
      names(inner) <- names(x)[i]
    }

    result <- c(result, inner)

  }

  result


}

