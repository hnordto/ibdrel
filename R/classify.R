prepareFeatures = function(lengthData, cutoff = 0) {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

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

classify = function(obs, pdfuns) {
  logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, log = T))
  res = exp(logprobs - matrixStats::logSumExp(logprobs))
  sort(res, decreasing = T)
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
