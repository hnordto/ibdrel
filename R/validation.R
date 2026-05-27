
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
