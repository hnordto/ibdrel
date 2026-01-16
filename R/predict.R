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


  return (pred)
}

# Helper function to retrieve each observation when multiple are inputted
get_inner_element <- function(x) {
  if (!is.list(x)) {
    return(list(x))
  } else {
    return(unlist(lapply(x, get_inner_element), recursive = FALSE))
  }
}
