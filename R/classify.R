prepareFeatures = function(lengthData, cutoff = 0) {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  # Remove empty lists (simulations with no segments > cutoff)
  lengthData <- lengthData[lapply(lengthData, length) > 0]

  features <- list(countPdf = lengths(lengthData),
                   totalPdf = vapply(lengthData, sum, FUN.VALUE = 1),
                   medianPdf = vapply(lengthData, median, FUN.VALUE = 1),
                   longestPdf = vapply(lengthData, max, FUN.VALUE = 1),
                   shortestPdf = vapply(lengthData, min, FUN.VALUE = 1))

  return (features)

}

preparePdfs = function(lengthData, cutoff = 0) {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  # Remove empty lists (simulations with no segments > cutoff)
  lengthData <- lengthData[lapply(lengthData, length) > 0]

  pdfs <- list(countPdf = lengths(lengthData) |>
                       density(from = 0) |> approxfun(rule = 2),
                     totalPdf = vapply(lengthData, sum, FUN.VALUE = 1) |>
                       density(from = 0) |> approxfun(rule = 2),
                     lengthPdf = unlist(lengthData, use.names = FALSE) |>
                       density(from = 0) |> approxfun(rule = 2),
                     medianPdf = vapply(lengthData, median, FUN.VALUE = 1) |>
                       density(from = 0) |> approxfun(rule = 2),
                     longestPdf = vapply(lengthData, max, FUN.VALUE = 1) |>
                       density(from = 0) |> approxfun(rule = 2),
                     shortestPdf = vapply(lengthData, min, FUN.VALUE = 1) |>
                       density(from = 0) |> approxfun(rule = 2))

  return (pdfs)

}

classProb = function(obs, pdfs, log = T) {
  p0 = pdfs$countPdf(length(obs))
  p1 = pdfs$totalPdf(sum(obs))
  p2 = pdfs$medianPdf(median(obs))
  p3 = pdfs$longestPdf(max(obs))
  p4 = pdfs$shortestPdf(min(obs))

  probs = c(p0,p1,p2,p3,p4)

  if (log) sum(log(probs)) else prod(probs)

}

classProbs = function(obs, pdfuns) {
  logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, log = T))
  logprobs
}

normalizeClassProbs = function(logprobs) {
  res = exp(logprobs - matrixStats::logSumExp(logprobs))
}

classify = function(obs, pdfuns, sort = TRUE) {
  logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, log = T))
  res = exp(logprobs - matrixStats::logSumExp(logprobs))

  if(isTRUE(sort)) {
    sort(res, decreasing = T)
  } else {
    res
  }
}


# Goodness-of-fit --------------------------

computeCovariance <- function(features) {
  cov(as.data.frame(features), use = "pairwise.complete.obs")
}

obsToFeatures <- function(obs) {
  obs.lst = list()
  obs.lst$obs[[1]] = obs # Same list structure as training data

  obs.features = lapply(obs.lst, prepareFeatures)

  obs.features
}

computeMahalanobis <- function(features, obs) {

  features <- as.data.frame(features)
  obs.features <- as.data.frame(obsToFeatures(obs))

  if (nrow(unique(features)) == 1) { # In the case of lineal of degree 1
    return (NA)
  }

  pca <- prcomp(features, center = TRUE, scale. = TRUE)
  k <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.95)[1]

  features.reduced <- pca$x[,1:k]

  obs.features.pca <- scale(obs.features,
                            center = pca$center,
                            scale = pca$scale) %*% pca$rotation[,1:k]

  covmat <- computeCovariance(features.reduced)

  if (all(covmat) == 0) { # If PO, covariance is 0 (no variation, always 50% IBD)
    return (NA)
  }

  features.colmeans <- colMeans(features.reduced)

  mdist <- stats::mahalanobis(obs.features.pca, features.colmeans, covmat, tol = 1e-20)

  mdist

}

distance <- function(obs, features) {
  res = unlist(sapply(features, function(x) computeMahalanobis(x, obs)))
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
    ped.rel = names(segmentDataTest)[i]
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
      prediction_top = classify(segment, pdfs, sort = TRUE)[1]

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
    colnames(res_mat) <- class_lst
    res <- res_mat
  }

  return (res)
} # Can be used to test both aggregation before and aggregation after!


trueClasses <- function(testsegments,
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
