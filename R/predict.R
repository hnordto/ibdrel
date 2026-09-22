#'
#'@export
predict <- function(observed, model,
                    normalize = TRUE, sort = TRUE) {

  pdfs = model$pdfs

  obs = get_inner_element(observed)

  if (length(obs) > 1 && isTRUE(sort)) {
    stop("Cannot sort multiple observations, as ranks may differ. Set sort = FALSE.")
  }

  pred = lapply(pdfs, function(pdf) {
    res_lst = lapply(obs, function(obs_vec) {
      classify(obs_vec, pdf, cutoff = model$cutoff, sort = sort)
    })

    do.call(rbind, res_lst)
  })


  if (normalize) {
    pred = lapply(pred, function(mat) {
      t(apply(mat, 1, normalizeLikelihoods))
    })
  }

  pred <- lapply(RESOLUTIONS, function(x) {
    aggregateProbs(pred, pedsMetadata(model$peds), x, TRUE)
  })
  names(pred) <- RESOLUTIONS


  return (pred)
}

aggregateProbs <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {
  lookup = lookupClass(NULL, metadata.agg.column, "rel", metadata, collapseDistant)

  if (isTRUE(collapseDistant)) {
    lookup["distant"] = "distant"
  }

  agg.ids = lookup[colnames(pred$eqclass.detailed)]

  t(rowsum(t(pred$eqclass.detailed), group = agg.ids))
}

aggregateProbs2 <- function(pred, metadata, metadata.agg.column, collapseDistant = TRUE) {
  lookup = lookupClass(NULL, metadata.agg.column, "eqclass.detailed", metadata, collapseDistant)

  group <- factor(lookup[colnames(pred$eqclass.detailed)], levels = unique(lookup[colnames(pred$eqclass.detailed)]))

  M <- matrix(0, nrow = ncol(pred$eqclass.detailed), ncol = nlevels(group), dimnames = list(colnames(pred$eqclass.detailed), levels(group)))

  M[cbind(seq_along(group), as.integer(group))] <- 1

  pred$eqclass.detailed %*% M
}

# Helper function to retrieve each observation when multiple are inputted
get_inner_element <- function(x) {
  if (!is.list(x)) {
    return(list(x))
  } else {
    return(unlist(lapply(x, get_inner_element), recursive = FALSE))
  }
}
