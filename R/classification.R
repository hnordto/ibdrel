#'
#'@export
prepareFeatures = function(lengthData, featureSel = c("count", "total"), cutoff = 0) {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  features <- list(
    count = if ("count" %in% featureSel) lengths(lengthData) else NULL,
    total = if ("total" %in% featureSel) vapply(lengthData, sum, FUN.VALUE = 1) else NULL,
    median = if ("median" %in% featureSel) vapply(lengthData, safe_median, FUN.VALUE = 1),
    longest = if ("longest" %in% featureSel) vapply(lengthData, safe_max, FUN.VALUE = 1) else NULL,
    shortest = if ("shortest" %in% featureSel) vapply(lengthData, safe_min, FUN.VALUE = 1) else NULL
  )

  features <- Filter(Negate(is.null), features) # Remove NULL (not used features)

  return (features)

}

# Helpers

safe_max <- function(x, default = 0) {
  if (length(x) == 0) default else max(x)
}

safe_min <- function(x, default = 0) {
  if (length(x) == 0) default else min(x)
}

safe_lengths <- function(lst) {
  sapply(lst, function(x) {
    if (identical(x,0)) 0 else length(x)
  })
}


#'
#'@export
preparePdfs = function(lengthData, featureSel = c("count", "total"), cutoff = 0, bw = "nrd0") {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  pdfs <- list(
    count = if ("count" %in% featureSel) lengths(lengthData) |>
      density(from = 0, bw = "nrd0", kernel = "rectangular") |> approxfun(rule = 2) else NULL,
    total = if ("total" %in% featureSel) vapply(lengthData, sum, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    length = if ("length" %in% featureSel) unlist(lengthData, use.names = FALSE) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    longest = if ("longest" %in% featureSel) longestPdf = vapply(lengthData, safe_max, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    shortest = if ("shortest" %in% featureSel) vapply(lengthData, safe_min, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL
  )

  pdfs <- Filter(Negate(is.null), pdfs)

  return (pdfs)

}

classProb = function(obs, pdfs, obsType = "segments", log = T) {

  var.names <- names(pdfs)

  probs = c()

  if (obsType == "segments") {
    if ("count" %in% var.names) {
      probs = c(probs, pdfs$count(length(obs)))
    }

    if ("total" %in% var.names) {
      probs = c(probs, pdfs$total(sum(obs)))
    }

    if ("length" %in% var.names) {
      probs = c(probs, pdfs$length(obs))
    }

    if ("longest" %in% var.names) {
      probs = c(probs, pdfs$longest(safe_max(obs)))
    }

    if ("shortest" %in% var.names) {
      probs = c(probs, pdfs$shortest(safe_min(obs)))
    }
  } else if (obsType == "summaries") {
    if ("count" %in% var.names) {
      probs = c(probs, pdfs$count(obs$count))
    }

    if ("total" %in% var.names) {
      probs = c(probs, pdfs$total(obs$total))
    }

    if ("longest" %in% var.names) {
      probs = c(probs, pdfs$longest(obs$max))
    }

    if ("shortest" %in% var.names) {
      probs = c(probs, pdfs$shortest(obs$min))
    }
  } else {
    stop("Invalid obsType. Must be either 'segments' (default) or 'summaries'.")
  }


#  p0 = pdfs$countPdf(length(obs))
#  p1 = pdfs$totalPdf(sum(obs))
#  p2 = pdfs$medianPdf(median(obs))
#  p3 = pdfs$longestPdf(max(obs))
#  p4 = pdfs$shortestPdf(min(obs))
  # p5 = pdfs$lengthPdf(obs)

 # probs = c(p0,p1,p2,p3,p4, p5)

  if (log) sum(log(probs)) else prod(probs)

}

classProbs = function(obs, pdfuns) {
  logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, log = T))
  logprobs
}

#'
#'@export
classify = function(obs, pdfuns, obsType = "segments", cutoff = 0, sort = TRUE) {

  if(obsType == "segments") {
    obs = obs[obs >= cutoff]

    if (identical(obs, numeric(0)) || is.null(obs)) {
      logprobs = rep(-Inf, length(pdfuns)) # Likelihood = 0 for all classes if obs < cutoff
      #} else if (sum(obs)/3391 < 1/2048) {
      #logprobs = rep(-Inf, length(pdfuns))
    } else {
      logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, obsType = "segments", log = T))
    }
  } else if (obsType == "summaries") {
    logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, obsType = "summaries", log = T))
  } else {
    stop("Invalid obsType. Must be either 'segments' (default) or 'summaries'.")
  }





  if(isTRUE(sort)) {
    sort(logprobs, decreasing = T)
  } else {
    logprobs
  }
}

#'
#' @export
normalizeLikelihoods = function(likelihoods) {
  exp(likelihoods - matrixStats::logSumExp(likelihoods))
}

# Goodness-of-fit --------------------------



computeCovariance <- function(features) {
  cov(as.data.frame(features), use = "pairwise.complete.obs")
}

#'
#'@export
obsToFeatures <- function(obs, cutoff, featureSel, isFeatures = FALSE) {

  if (isFALSE(isFeatures)) {
    obs.lst = list()
    obs.lst$obs[[1]] = obs # Same list structure as training data

    obs.features = lapply(obs.lst, prepareFeatures, cutoff = cutoff, featureSel = featureSel)
  } else {
    obs.features = list("obs" = obs)
  }

  obs.features
}

computeMahalanobis <- function(features, obs, cutoff, featureSel, threshold, isFeatures) {

  features <- as.data.frame(features)
  obs.features <- as.data.frame(obsToFeatures(obs, cutoff, featureSel, isFeatures))

  if (nrow(unique(features)) == 1) { # In the case of lineal of degree 1
    return (list(dist = NA, threshold = NA))
  }

  #pca <- prcomp(features, center = TRUE, scale. = TRUE)
  #k <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.95)[1]

  #features.reduced <- pca$x[,1:k]

  #obs.features.pca <- scale(obs.features,
  #                          center = pca$center,
  #                          scale = pca$scale) %*% pca$rotation[,1:k]

  covmat <- computeCovariance(features)

  if (all(covmat == 0)) { # If PO, covariance is 0 (no variation, always 50% IBD)
    return (list(dist = NA, threshold = NA))
  }

  features.colmeans <- colMeans(features)

  colnames(obs.features) <- colnames(features)

  mdist <- stats::mahalanobis(rbind(features, obs.features), features.colmeans, covmat, tol = 1e-20)

  list(dist = tail(mdist, 1), threshold = quantile(mdist, threshold, names = FALSE))

}

#'
#'@export
distance <- function(obs, features, cutoff, featureSel, threshold, isFeatures) {

  #order.idx <- match(names(orderLst), names(features))

  #features <- features[order.idx]

  res = lapply(features, function(x) computeMahalanobis(x, obs, cutoff, featureSel, threshold, isFeatures))

  res
}

#'
#'@export
trueClasses <- function(testsegments,
                        metadata,
                        agg.level) {

  truth <- c()

  for (i in 1:length(testsegments)) {
    true = metadata |>
      dplyr::filter(rel == names(testsegments)[i]) |>
      dplyr::select(agg.level) |>
      as.character()

    truth <- c(truth, rep(true, length(testsegments[[i]])))
  }

  return (truth)
}
