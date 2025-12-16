
predict <- function(observed, model,
                    normalize = TRUE, sort = TRUE) {

  pdfs = model$pdfs

  pred = lapply(pdfs, function(x) {
    classify(observed, x, cutoff = model$cutoff, sort = sort)
  })

  if (normalize) {
    pred = lapply(pred, function(x) {
      normalizeLikelihoods(x)
    })
  }

  return (pred)
}

