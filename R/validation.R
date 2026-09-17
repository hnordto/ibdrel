
topPred = function(pred, n = NULL) {
  if (is.null(n)) {
    n = nrow(pred[[1]])
  }

  tp = lapply(pred, function(mat) {
    (apply(mat, 1, function(row) {
      colnames(mat)[order(row, decreasing = TRUE)[1:n]]
    }))
  })

  return (tp)
}

#'
#'@export
sensitivity <- function(pred, true, k=1) {
  classes <- colnames(pred)
  true <- factor(true, levels = classes)
  true_chr <- as.character(true)

  res <- lapply(k, function(kk) {

    topk_idx <- t(apply(pred, 1, function(x) {
      order(x, decreasing = TRUE)[1:kk]
    }))

    topk_class <- matrix(
      classes[topk_idx],
      nrow = nrow(pred),
      ncol = kk
    )

    res.k <- lapply(classes, function(c) {
      is_true_c <- (true_chr == c)

      is_pred_in_top_k <- apply(topk_class, 1, function(x) {
        c %in% x
      })

      TP <- sum(is_true_c & is_pred_in_top_k)
      FN <- sum(is_true_c & !is_pred_in_top_k)

      sens <- if ((TP+FN) > 0) TP / (TP+FN) else NA_real_

      data.frame(
        class = c,
        sensitivity = sens,
        k = kk
      )

    })

    do.call(rbind, res.k)

  })

  do.call(rbind, res)

}

ppv <- function(pred, true) {
  classes <- colnames(pred)
  true <- factor(true, levels = classes)
  true_chr <- as.character(true)


}

model.logloss <- function(segments.test, features, cutoff) {
  metadata = pedsMetadata(ibdrel_unilineal$peds)

  true.eqdetailed <- factor(trueClasses(segments.test, metadata, "eqclass.detailed"))
  true.eq <- factor(trueClasses(segments.test, metadata, "eqclass"))
  true.deg <- factor(trueClasses(segments.test, metadata, "degree"))

  res <- data.frame(features = character(),
                    logloss.eqdetailed = numeric(),
                    logloss.eq = numeric(),
                    logloss.eeg = numeric())

  for (i in 1:length(features)) {
    feature <- features[[i]]

    model = fitModel(features = feature, cutoff = cutoff)

    pred = predict(segments.test, model, sort = FALSE)

    pred.eqdetailed <- pred$eqclass.detailed
    pred.eq <- pred$eqclass
    pred.deg <-pred$deg

    true.eqdetailed <- factor(true.eqdetailed, levels = colnames(pred.eqdetailed))
    true.eq <- factor(true.eq, levels = colnames(pred.eq))
    true.deg <- factor(true.deg, levels = colnames(pred.deg))

    w.eqdetailed <- 1 / table(true.eqdetailed)[as.character(true.eqdetailed)]
    w.eq <- 1 / table(true.eq)[as.character(true.eq)]
    w.deg <- 1 / table(true.deg)[as.character(true.deg)]

    logloss.eqdetailed <- SLmetrics::weighted.logloss(actual = true.eqdetailed, response = pred.eqdetailed, w = w.eqdetailed)
    logloss.eq <- SLmetrics::weighted.logloss(actual = true.eq, response = pred.eq, w = w.eq)
    logloss.deg <- SLmetrics::weighted.logloss(actual = true.deg, response = pred.deg, w = w.deg)

    res <- rbind(res, data.frame(features = toString(feature),
                                 logloss.eqdetailed = logloss.eqdetailed,
                                 logloss.eq = logloss.eq,
                                 logloss.deg = logloss.deg))
  }

  return (res)
}
