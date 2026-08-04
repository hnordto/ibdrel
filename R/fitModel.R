#' Fit a relationship classification model
#'
#' Train a non-parametric Naive Bayes classification model on IBD segment lengths.
#'
#' A probability density function (PDF) for each combination of class and feature
#' is estimated by kernel density estimation (KDE). Segments below the specified
#' threshold (`cutoff`) are treated as no IBD.
#'
#' @param data Nested list of training data.
#'   A hierarchical list structured as:
#'   \code{data[[class_level]][[class]] -> list of training instances},
#'   where:
#'   \itemize{
#'     \item \strong{class_level}: grouping of relationship classes
#'     \item \strong{class}: a specific relationship
#'     \item each element: a training instance (feature vector)
#'   }
#'   Allowed class levels and classes depend on your dataset and should match the expected hierarchy.
#'
#' @param features Character vector specifying which features to use.
#'   Allowed values: \code{"count"}, \code{"total"}, \code{"median"}, \code{"longest"}, \code{"shortest"}.
#'
#' @param cutoff Numeric value indicating the minimum segment length to be considered IBD.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \strong{pdfs}: nested list of PDFs, following the structure of the input \code{data}
#'     \item \strong{cutoff}: the cutoff used for IBD segment filtering
#'     \item \strong{features}: the features used for the model
#'   }
#'
#' @examples
#' ### Example 1: Fit default model ###
#' model <- fitModel()
#' model
#'
#' @export
fitModel <- function(data = NULL,
                     features = "default",
                     cutoff = 7,
                     collapseDistant = TRUE) {
  if (is.null(data)) {
    data = loadData(collapseDistant = collapseDistant)
  }

  if (all(features == "default")) {
    pdfs = lapply(names(data), function(nm) {
      x <- data[[nm]]
      lapply(x, preparePdfs, featureSel = featureSelection.default[[nm]],
             cutoff = cutoff)
    })
    names(pdfs) <- names(data)
  } else {
    if (any(!(features %in% VALID.FEATURES))) {
      stop("Invalid feature(s). Must be one or multiple of ",
           paste(VALID.FEATURES, collapse = ", "))
    }

    pdfs = lapply(data, function(x) {
      lapply(x, preparePdfs,
             featureSel = features,
             cutoff = cutoff)
    })

  }

  model = list(
    pdfs = pdfs,
    cutoff = cutoff,
    features = features
  )

  return (model)

}

featureSelection.default = list("eqclass.detailed" = c("count", "length"),
                                "eqclass" = c("count", "total"),
                                "kappa" = c("count", "total"),
                                "kinship" = c("count", "total"),
                                "degree" = c("count", "total"))
