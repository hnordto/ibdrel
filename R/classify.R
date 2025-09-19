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

safe_median <- function(x, default = 0) {
  if(length(x) == 0) default else median(x)
}

safe_lengths <- function(lst) {
  sapply(lst, function(x) {
    if (identical(x,0)) 0 else length(x)
  })
}

preparePdfs = function(lengthData, featureSel = c("count", "total"), cutoff = 0, bw = "nrd0") {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  pdfs <- list(
    count = if ("count" %in% featureSel) lengths(lengthData) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    total = if ("total" %in% featureSel) vapply(lengthData, sum, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    median = if ("median" %in% featureSel) vapply(lengthData, safe_median, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    longest = if ("longest" %in% featureSel) longestPdf = vapply(lengthData, safe_max, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL,
    shortest = if ("shortest" %in% featureSel) vapply(lengthData, safe_min, FUN.VALUE = 1) |>
      density(from = 0, bw = bw) |> approxfun(rule = 2) else NULL
  )

  pdfs <- Filter(Negate(is.null), pdfs)

  return (pdfs)

}

classProb = function(obs, pdfs, log = T) {

  var.names <- names(pdfs)

  probs = c()

  if ("count" %in% var.names) {
    probs = c(probs, pdfs$count(length(obs)))
  }

  if ("total" %in% var.names) {
    probs = c(probs, pdfs$total(sum(obs)))
  }

  if ("median" %in% var.names) {
    probs = c(probs, pdfs$median(safe_median(obs)))
  }

  if ("longest" %in% var.names) {
    probs = c(probs, pdfs$longest(safe_max(obs)))
  }

  if ("shortest" %in% var.names) {
    probs = c(probs, pdfs$shortest(safe_min(obs)))
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

normalizeClassProbs = function(logprobs) {
  res = exp(logprobs - matrixStats::logSumExp(logprobs))
}

classify = function(obs, pdfuns, cutoff = 0, sort = TRUE) {
  obs = obs[obs >= cutoff]

  if (identical(obs, numeric(0)) || is.null(obs)) {
    logprobs = rep(-Inf, length(pdfuns)) # Likelihood = 0 for all classes if obs < cutoff
  } else {
    logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, log = T))
  }



  if(isTRUE(sort)) {
    sort(logprobs, decreasing = T)
  } else {
    logprobs
  }
}

normalizeLikelihoods = function(likelihoods) {
  exp(likelihoods - matrixStats::logSumExp(likelihoods))
}

# Goodness-of-fit --------------------------



computeCovariance <- function(features) {
  cov(as.data.frame(features), use = "pairwise.complete.obs")
}

obsToFeatures <- function(obs, cutoff, featureSel) {
  obs.lst = list()
  obs.lst$obs[[1]] = obs # Same list structure as training data

  obs.features = lapply(obs.lst, prepareFeatures, cutoff = cutoff, featureSel = featureSel)

  obs.features
}

computeMahalanobis <- function(features, obs, cutoff, featureSel) {

  features <- as.data.frame(features)
  obs.features <- as.data.frame(obsToFeatures(obs, cutoff, featureSel))

  if (nrow(unique(features)) == 1) { # In the case of lineal of degree 1
    return (NA)
  }

  #pca <- prcomp(features, center = TRUE, scale. = TRUE)
  #k <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.95)[1]

  #features.reduced <- pca$x[,1:k]

  #obs.features.pca <- scale(obs.features,
  #                          center = pca$center,
  #                          scale = pca$scale) %*% pca$rotation[,1:k]

  covmat <- computeCovariance(features)

  if (all(covmat == 0)) { # If PO, covariance is 0 (no variation, always 50% IBD)
    return (NA)
  }

  features.colmeans <- colMeans(features)

  mdist <- stats::mahalanobis(obs.features, features.colmeans, covmat, tol = 1e-20)

  mdist

}

computeLOF <- function(features, obs) {
  features <- as.data.frame(features)
  obs.features <- as.data.frame(obsToFeatures(obs))
  colnames(obs.features) <- colnames(features)

  data <- rbind(features, obs.features)
  #data$countPdf <- log1p(data$countPdf)
  data <- scale(data)

  # Use a subsample if nrow > 1000 (approximation to enhance efficiency)

  if (nrow(data) > 10000) {
    idx <- sample(1:nrow(data), 10000, replace = FALSE)
    data <- data[idx,]
  }

  lofs <- dbscan::lof(data, minPts = 30)
  lof.threshold <- quantile(lofs, p = .9)
  lof <- lofs[length(lofs)] # LOF of obs

  list(lof = lof, threshold = lof.threshold)



}

LOF <- function(obs, features, orderLst, top_n) {

  order.idx <- match(names(orderLst), names(features))

  features.full <- features
  features <- features[order.idx]
  features <- features[1:top_n]

  res = lapply(features, function(x) computeLOF(x, obs))

  lof <- as.numeric(lapply(res, function(r) r$lof))
  threshold <- as.numeric(lapply(res, function(r) r$threshold))

  lof <- c(lof, rep(NA, length(features.full)-top_n))
  threshold <- c(threshold, rep(NA, length(features.full)-top_n))

  list(lof = lof, threshold = threshold)

}

distance <- function(obs, features, cutoff, featureSel) {

  #order.idx <- match(names(orderLst), names(features))

  #features <- features[order.idx]

  res = unlist(sapply(features, function(x) computeMahalanobis(x, obs, cutoff, featureSel)))

  res
}



## Performance

testClassifier <- function(testsegments,
                           pdfs,
                           metadata,
                           agg.level,
                           all = FALSE) {
  if(isFALSE(all)) {
    res <- data.frame(true = character(),
                      pred = character(),
                      prob = numeric())
  } else {
    res_mat <- matrix(data = NA,
                      nrow = sum(sapply(testsegments, length)),
                      ncol = nrow(unique(metadata[agg.level])))
  }

  class_lst = c()
  k = 1

  for (i in 1:length(testsegments)) {
    ped.rel = names(testsegments)[i]
    segments = testsegments[[i]]

    true = metadata |>
      filter(rel == ped.rel) |>
      select(agg.level) |>
      as.character()

    if (!(true %in% class_lst)) {
      class_lst = c(class_lst, true)
    }


    for (segment in segments) {
      prediction = classify(segment, pdfs, sort = FALSE) # Very important!
      prediction = normalizeClassProbs(prediction)
      prediction_top = classify(segment, pdfs, sort = TRUE)[1]
    #  prediction_top = normalizeClassProbs(prediction_top)

      pred = metadata |>
        filter(rel == names(prediction_top)) |>
        select(agg.level) |>
        as.character()

      if (isFALSE(all)) {
        res.tmp <- data.frame(true = true,
                              pred = pred,
                              prob = prediction_top)
        res <- rbind(res, res.tmp)
      } else {
        res_mat[k, ] = prediction
        k = k +1
      }
    }
  }

  if (isFALSE(all)) {
    res <- res
  } else {
    colnames(res_mat) <- names(prediction)
    res <- res_mat
  }

  return (res)
} # Can be used to test both aggregation before and aggregation after!


trueClasses <- function(testsegments,
                        metadata,
                        agg.level) {

  truth <- c()

  for (i in 1:length(testsegments)) {
    true = metadata |>
      filter(rel == names(testsegments)[i]) |>
      select(agg.level) |>
      as.character()

    truth <- c(truth, rep(true, length(testsegments[[i]])))
  }

  return (truth)
}
