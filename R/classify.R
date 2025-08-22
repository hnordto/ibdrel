prepareFeatures = function(lengthData, cutoff = 0) {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  # Remove empty lists (simulations with no segments > cutoff)
  #lengthData <- lengthData[lapply(lengthData, length) > 0]

  features <- list(countPdf = lengths(lengthData),
                   totalPdf = vapply(lengthData, sum, FUN.VALUE = 1))
                   #medianPdf = vapply(lengthData, safe_median, FUN.VALUE = 1),
                   #longestPdf = vapply(lengthData, safe_max, FUN.VALUE = 1),
                   #shortestPdf = vapply(lengthData, safe_min, FUN.VALUE = 1))

  return (features)

}

# Helpers

safe_max <- function(x, default = 0) {
  val <- max(x)
  if (length(x) == 0 || val == -Inf) default else val
}

safe_min <- function(x, default = 0) {
  val <- min(x)
  if (length(x) == 0 || val == Inf) default else val
}

safe_median <- function(x, default = 0) {
  val <- median(x)
  if(length(x) == 0 || is.na(val)) default else val
}

safe_lengths <- function(lst) {
  sapply(lst, function(x) {
    if (identical(x,0)) 0 else length(x)
  })
}

preparePdfs = function(lengthData, cutoff = 0, bw = "nrd0") {
  if (cutoff > 0) {
    lengthData = lapply(lengthData, function(v) v[v >= cutoff])
  }

  # Remove empty lists (simulations with no segments > cutoff)
  #lengthData <- lengthData[lapply(lengthData, length) > 0]

  pdfs <- list(countPdf = lengths(lengthData) |>
                       density(from = 0, bw = bw) |> approxfun(rule = 2),
                     totalPdf = vapply(lengthData, sum, FUN.VALUE = 1) |>
                       density(from = 0, bw = bw) |> approxfun(rule = 2),
                     #lengthPdf = unlist(lengthData, use.names = FALSE) |>
                    #   density(from = 0) |> approxfun(rule = 2),
                     medianPdf = vapply(lengthData, safe_median, FUN.VALUE = 1) |>
                       density(from = 0, bw = bw) |> approxfun(rule = 2),
                     longestPdf = vapply(lengthData, safe_max, FUN.VALUE = 1) |>
                       density(from = 0, bw = bw) |> approxfun(rule = 2),
                     shortestPdf = vapply(lengthData, safe_min, FUN.VALUE = 1) |>
                       density(from = 0, bw = bw) |> approxfun(rule = 2))

  return (pdfs)

}

classProb = function(obs, pdfs, log = T) {

  var.names <- names(pdfs)

  probs = c()

  if ("countPdf" %in% var.names) {
    probs = c(probs, pdfs$countPdf(length(obs)))
  }

  if ("totalPdf" %in% var.names) {
    probs = c(probs, pdfs$totalPdf(sum(obs)))
  }

  if ("lengthPdf" %in% var.names) {
    probs = c(probs, pdfs$lengthPdf(obs))
  }

  if ("medianPdf" %in% var.names) {
    probs = c(probs, pdfs$medianPdf(median(obs)))
  }

  if ("longestPdf" %in% var.names) {
    probs = c(probs, pdfs$longestPdf(max(obs)))
  }

  if ("shortestPdf" %in% var.names) {
    probs = c(probs, pdfs$shortestPdf(min(obs)))
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
  obs = obs[obs > cutoff]
  logprobs = sapply(pdfuns, function(pdfs) classProb(obs, pdfs, log = T))

  if(isTRUE(sort)) {
    sort(logprobs, decreasing = T)
  } else {
    logprobs
  }
}

normalizePosteriors = function(posteriors) {
  exp(posteriors - matrixStats::logSumExp(posteriors))
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

  #pca <- prcomp(features, center = TRUE, scale. = TRUE)
  #k <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.95)[1]

  #features.reduced <- pca$x[,1:k]

  #obs.features.pca <- scale(obs.features,
  #                          center = pca$center,
  #                          scale = pca$scale) %*% pca$rotation[,1:k]

  covmat <- computeCovariance(features)

  if (all(covmat) == 0) { # If PO, covariance is 0 (no variation, always 50% IBD)
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

distance <- function(obs, features, orderLst, top_n) {

  order.idx <- match(names(orderLst), names(features))

  features.full <- features
  features <- features[order.idx]
  features <- features[1:top_n]

  res = unlist(sapply(features, function(x) computeMahalanobis(x, obs)))

  res = c(res, rep(NA, length(features.full)-top_n))

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
