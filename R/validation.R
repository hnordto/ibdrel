
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
